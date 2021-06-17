using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO,
      StatsBase,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi


###
### Gurobi licence fidgeting (have it do only one licence check)
###

if !(@isdefined env)
    const env = Gurobi.Env()
end


###
### Define graph and network
###

# Load simple graph
g = loadgraph("g_test.mg", MGFormat())
positions = get_node_pos(g)
# Make directed
g_directed = DiGraph(g.graph)
sorted_edges = collect(edges(g_directed))

# Make road network
a = edge_lengths(g_directed, positions)
b = resource_allocation(g_directed, a)
edge_params = Dict(:a => a, :b => b)
node_params = Dict(:pos => positions)
rn = RoadNetwork(g_directed, edge_params, node_params)

n = nv(rn.g)
m = ne(rn.g)

# Demand structure
ods = [(1,9)]
n_ods = length(ods)
demands = [1.0]
γ = 0.1

###
### Find minimum spanning tree
###

mst = kruskal_mst(g)
mst_rev = reverse.(mst)

mst_edges = [(e.src, e.dst) for e in mst]
mst_rev_edges = [(e.src, e.dst) for e in mst_rev]

# Identify MST edges
mst_indices = findfirst.(isequal.(mst), [sorted_edges])
mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])

protected_indices = vcat(mst_indices, mst_rev_indices)
protected_edges = vcat(mst, mst_rev)

# Edges not in MST
other_indices = [] 
other_edges = []
for (i, ed) in enumerate(sorted_edges)
    if !in(ed, protected_edges)
        push!(other_indices, i)
        push!(other_edges, ed)
    end     
end

#
# MST graph (most convenient way of doing things with it)
#

g_mst = DiGraph(protected_edges)

# Dist matrix on MST graph
mst_e_lens = rn.edge_params[:a][sort(protected_indices)]
mst_distmx = spzeros(n,n) # Annoying way of using LightGraphs functions
for (k,ed) in enumerate(edges(g_mst)) 
    i, j = ed.src, ed.dst
    mst_distmx[i,j] = mst_e_lens[k]      
end
ds = dijkstra_shortest_paths(g_mst, ods[1][1], mst_distmx)
tree_path = enumerate_paths(ds, ods[1][2])




### 
### Save for plotting
###

#save_graph_tikz_edg(g_directed,
#                    [mst_edges, mst_rev_edges],
#                    positions,
#                    "test_net_mst.tex",
#                    imp_nodes=[ods[1][1], ods[1][2]],
#                    imp_labels=["O", "D"])




###
### Choice of restricted roads
###

# Start with random choice
max_res = m - length(mst)
# Random choice (proportion)
prop = 0.1
n_restricted = Int(floor(prop*max_res))
#bounded_indices = sample(other_indices, n_restricted)
bounded_indices = other_indices

###
### Define Traffic Assignment
###

# Cost structure
A = my_incidence_matrix(rn.g)
a = rn.edge_params[:a]
B = diagm(rn.edge_params[:b])

# Vectors for conservation constraints
d_vects_hv = SparseVector{Float64,Int64}[]
d_vects_av = SparseVector{Float64,Int64}[]
for i in 1:length(ods)
    s, d = ods[i]
    d_vec = spzeros(n)
    d_vec[s] = -demands[i]
    d_vec[d] = demands[i]
    
    push!(d_vects_hv, (1-γ).*d_vec)
    push!(d_vects_av, γ.*d_vec)
end

# Objective and penalty
me_obj_func(x, y) = sum( B[i,i]*(x[i]+y[i])^2 + 
                             a[i]*x[i] + 0.5*a[i]*y[i] for i in 1:length(x))

penalty_function(x, r, m) = r*(x^m)
function bound_penalty(x, penalty_param, indices)
    suma = 0
    for i in indices
        suma += penalty_function(x[i], 0.5*a[i], 2)
    end
    penalty_param * suma 
end                              

penalised_objective(x, y, penalty_param, bound_indices) = me_obj_func(x, y) + bound_penalty(x, penalty_param, bound_indices)

###
### JuMP model construction
###

"""
For now, this function just wraps sequence of problems 
for the penalty method to make interactive use easier
"""
function penalty_stap_wrapper(r_array)

    #stap = Model(Ipopt.Optimizer)
    #stap = Model(Gurobi.Optimizer)
    stap = Model(() -> Gurobi.Optimizer(env))

    # OD specific link flows HV         
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)

    # OD specific link flows AV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    #hv_flow = sum(x, dims=2)
    #av_flow = sum(y, dims=2)

    # Conservation constraints
    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])


    # Unpenalised objective
    @objective(stap, Min, me_obj_func(x, y))
    
    ### 
    ### Solve unpenalised
    ###
    
    optimize!(stap)

    # Update
    sols_x = value.(x)
    sols_y = value.(y)

    iters_x = [sols_x[:]]
    iters_y = [sols_y[:]]
    
    ###
    ### Sequence of penalised sub-problems
    ###
    
    #r_array = [10, 100, 1000, 10000, 100000, 1000000]
    
    for r in r_array
        @objective(stap,
                   Min,
                   penalised_objective(x, y, r, bounded_indices))
        
        for k in 1:m
            set_start_value(x[k,1], sols_x[k])
            set_start_value(y[k,1], sols_y[k])
        end           
        
        optimize!(stap)
        
        sols_x = value.(x)
        sols_y = value.(y)

        push!(iters_x, sols_x[:])
        push!(iters_y, sols_y[:])
    end
    return iters_x, iters_y 
end

#r_array = [10, 100, 1000, 10000, 100000, 1000000]
r_array = [10^i for i in 1:10]

xx, yy = penalty_stap_wrapper(r_array)
agg = xx + yy

r_array2 = [10^5]

xx2, yy2 = penalty_stap_wrapper(r_array2)
agg2 = xx2 + yy2
