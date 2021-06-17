using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi

###
### Load simple graph
###

g = loadgraph("g_test.mg", MGFormat())
positions = get_node_pos(g)

g_directed = DiGraph(g.graph)
sorted_edges = collect(edges(g_directed))

# Make road network
a = edge_lengths(g_directed, positions)
b = resource_allocation(g_directed, a)
edge_params = Dict(:a => a, :b => b)
node_params = Dict(:pos => positions)
rn = RoadNetwork(g_directed, edge_params, node_params)

###
### Find minimum spanning tree
###

mst = kruskal_mst(g)
mst_rev = reverse.(mst)

mst_edges = [(e.src, e.dst) for e in mst]
mst_rev_edges = [(e.src, e.dst) for e in mst_rev]

# Edge indices (it might be better to work with the dictionaries)
mst_indices = findfirst.(isequal.(mst), [sorted_edges])
mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])

protected_indices = vcat(mst_indices, mst_rev_indices)
protected_edges = vcat(mst, mst_rev)

# Save graph with MST highlighted
#save_graph_tikz_edg(g, [mst_edges, mst_rev_edges], positions, "test_net_mst.tex")


###
### Define Traffic Assignment
###

ods = [(1,13)]
n_ods = length(ods)
demands = [1.0]
γ = 0.5

n = nv(rn.g)
m = ne(rn.g)

A = incidence_matrix(rn.g)
a = rn.edge_params[:a]
B = diagm(rn.edge_params[:b])

# demand vectors for conservation constraints.

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

me_obj_func(x, y, z) = sum( B[i,i]*(z[i]*x[i]+y[i])^2 + 
                             a[i]*z[i]*x[i] + y[i]*a[i]*0.5 for i in 1:length(x))

stap = Model(Gurobi.Optimizer)
set_silent(stap)

# OD specific link flows HV
@variable(stap, x[1:m,1:n_ods])
@constraint(stap, x .>= 0)

# OD specific link flows HV
@variable(stap, y[1:m,1:n_ods])
@constraint(stap, y .>= 0)

hv_flow = sum(x, dims=2)
av_flow = sum(y, dims=2)
#agg_flow = av_flow + hv_flow

@variable(stap, z[1:m, 1:n_ods], Bin)

for i in protected_indices
    @constraint(stap, z[i] == 1)
end

@constraint(stap,
        conservation_hvs[i in 1:n_ods],
        A*x[:,i] .== d_vects_hv[i])

@constraint(stap,
            conservation_avs[i in 1:n_ods],
            A*y[:,i] .== d_vects_av[i])



ex = @NLexpression(stap, sum( B[i,i]*(z[i]*x[i]+y[i])^2 + 
                             a[i]*z[i]*x[i] + y[i]*a[i]*0.5 for i in 1:length(x)))

@NLobjective(stap, Min, ex)
#@objective(stap,
#          Min,
#           mixed_objective) 

optimize!(stap)
