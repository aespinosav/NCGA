using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO,
      StatsBase,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi,
      Ipopt

# Fix incidence matrix in LightGraphs (for older versions of LightGraphs)
#import LightGraphs.LinAlg.incidence_matrix
function my_incidence_matrix(G)
    I = vcat([src(e) for e in edges(G)], [dst(e) for e in edges(G)])
    J = vcat(collect(1:ne(G)), collect(1:ne(G)))
    V = vcat(fill(-1, ne(G)), fill(1, ne(G)))
    return sparse(I, J, V)
end

###
### Gurobi licence check fidgeting
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
### Find MST
###


mst = kruskal_mst(g)
mst_rev = reverse.(mst)

mst_edges = [(e.src, e.dst) for e in mst]
mst_rev_edges = [(e.src, e.dst) for e in mst_rev]

# Identify MST edges (both ways)
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

# MST graph
g_mst = DiGraph(protected_edges)
# Dist matrix on MST graph
mst_e_lens = rn.edge_params[:a][sort(protected_indices)]
mst_distmx = spzeros(n,n) # Annoying way of using LightGraphs functions
for (k,ed) in enumerate(edges(g_mst)) 
    i, j = ed.src, ed.dst
    mst_distmx[i,j] = mst_e_lens[k]      
end
# shortest path on tree
ds = dijkstra_shortest_paths(g_mst, ods[1][1], mst_distmx)
tree_path = enumerate_paths(ds, ods[1][2])


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


function stap_wrapper(rn, ods; solver=:gurobi)

    # Structure
    n = nv(rn.g)
    m = ne(rn.g)
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

    # Objective (Mixed multi-class Beckmann UE/SO)
    me_obj_func(x, y) = sum( B[i,i]*(x[i]+y[i])^2 + 
                             a[i]*x[i] + 
                             0.5*a[i]*y[i] for i in 1:length(x))

    # Choose solver
    if solver == :gurobi
        stap = Model(() -> Gurobi.Optimizer(env))
    elseif solver == :ipopt
        stap = Model(Ipopt.Optimizer)
    end

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)

    # OD specific link flows AV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    # Conservation constraints
    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])

    #Bound constraints (HV flow set to 0 on some links)
    for i in bounded_indices
        for k in 1:n_ods
            @constraint(stap, x[i,k]==0)
        end
    end


    # Unpenalised objective
    @objective(stap, Min, me_obj_func(x, y))

    # Solve
    optimize!(stap)
    return value.(x)[:], value.(y)[:] 
end

xx, yy = stap_wrapper(rn, ods, solver=:ipopt)
agg = xx + yy

# Difference in total cost between using Gurobi and Ipopt is 3.9195797452151737e-7
