using TrafficNetworks2,
      SkeletonCities2,
      MHPaths,
      LightGraphs,
      MetaGraphs,
      SimpleWeightedGraphs,
      GraphIO,
      StatsBase,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi,
      Ipopt,
      Plots,
      DelimitedFiles

###
### Incidence matrix correction (and other functions)
###

include("tools.jl")
include("penalty.jl")

###
### Gurobi licence environment
###

if !(@isdefined env)
    const env = Gurobi.Env()
end

# network generated with net_gen_medium.jl (n=8, α_hat=0.8, β=1.5)
mg = loadgraph("g_test_medium.mg", MGFormat())
rn = skel2rn(mg)

a = rn.edge_params[:a]
b = rn.edge_params[:b]

sorted_edges = collect(edges(rn.g))

# STAP parameters

O = 1
D = 64

ods = [(O, D)]
d = 0.1

# ME parameters
γ_start = 0.0
γ_step = 0.1
γ_stop = 1.0

γ_array = γ_start:γ_step:γ_stop

n_ods = length(ods)
demands = [d]

# MH parameters
μ = 4
p_splice = 0.75

num_metropolis_steps = 5000
stagger = 5

# Penalty parameter array
r_array = [10^i for i in 1:7]


# Minimum Spanning Tree
mst = kruskal_mst(mg)
mst_rev = reverse.(mst)

mst_edges = [(e.src, e.dst) for e in mst]
mst_rev_edges = [(e.src, e.dst) for e in mst_rev]

# Identify MST edges
mst_indices = findfirst.(isequal.(mst), [sorted_edges])
mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])

# MST
pi_1 = vcat(mst_indices, mst_rev_indices)
pe_1 = vcat(mst, mst_rev)


# Metropolis-Hastings (Flötteröd + Bierlaire)
g_for_sampling = SimpleWeightedDiGraph(weighted_adjacency(mg))

mh = MHInstance(g_for_sampling, ods[1][1], ods[1][2], μ, p_splice)
mh_evolve!(mh, num_metropolis_steps, weight_func)

paths = unique([p.Γ for p in mh.history])
# stagger in steps still has to be properly selcted (this is a good "guess")
indep_paths = paths[1:stagger:end]

# paths are returned as lists of consecutive nodes
indep_edge_paths = Array{Edge,1}[]
for p in indep_paths
    ed_path = Array{Edge,1}(undef, length(p)-1)
    for (i, n) in enumerate(p[1:end-1])
        ed_path[i] = Edge(n, p[i+1])
    end
    push!(indep_edge_paths, ed_path)
end

indep_paths_ind = [findfirst.(isequal.(p), [sorted_edges])
                   for p in indep_edge_paths]

#writedlm("test_edge_indexed_paths.dat", indep_paths_ind)
#test_p = readdlm("test_edge_indexed_paths.dat")

# MH Paths
pi_2 = vcat(indep_paths_ind...)
pe_2 = vcat(indep_edge_paths...)

# HV edges (from MST and MH)
pi_both = unique(vcat(pi_1, pi_2))
pe_both = unique(vcat(pe_1, pe_2))

# AV exclusive edges
oi = []
oe = []
for (i, ed) in enumerate(sorted_edges)
    # MST
    if !in(ed, pe_both)
        push!(oi, i)
        push!(oe, ed)
    end
end

###
### Plot nets with paths and save paths in file
###

#save_paths_tikz(rn.g,
#                indep_paths,
#                get_node_pos(mg),
#                "tree_debug.tex",
#                imp_nodes=ods[1],
#                imp_labels=["O", "D"],
#                scale=10,
#                standalone_doc=true)
#
#edges_for_plotting = [(e.src, e.dst) for e in pe_both]
#save_graph_tikz_edg(rn.g,
#                    [edges_for_plotting],
#                    get_node_pos(mg),
#                    "tree_and_paths_debug.tex",
#                    edge_labels=false,
#                    imp_nodes=[1, 64],
#                    imp_labels=["O", "D"])
#
#writedlm("edge_ind_paths_debug.dat", indep_paths_ind)


###
### Assignment
###

### UE and SO

demand_array = 0.0001:0.001:0.3

UE_costs = zeros(length(demand_array))
SO_costs = zeros(length(demand_array))

for (k,d) in enumerate(demand_array)

    ue_flows = multi_pair_stap_nc(rn, ods, d, regime=:ue)
    so_flows = multi_pair_stap_nc(rn, ods, d, regime=:so)

    ue_cost = total_cost(ue_flows, a, b)
    so_cost = total_cost(so_flows, a, b)

    UE_costs[k] = ue_cost
    SO_costs[k] = so_cost
end

poa = UE_costs ./ SO_costs

#=
    For the example network at least, the PoA peak happens at d=0.101
    with three significant figures
=#

ue_flows = multi_pair_stap_nc(rn, ods, d, regime=:ue)
so_flows = multi_pair_stap_nc(rn, ods, d, regime=:so)

ue_cost = total_cost(ue_flows, a, b)
so_cost = total_cost(so_flows, a, b)

######################################################################


### ME with different penetration rates (and diff penalty method iters)
results = me_excluded_assignment_penalty_pr(rn,
                                            ods,
                                            demands,
                                            γ_array,
                                            oi,
                                            r_array)

results_hv, results_av, results_agg = results

end_results_hv = [results_hv[i][end] for i in 1:length(γ_array)]
end_results_av = [results_av[i][end] for i in 1:length(γ_array)]
end_results_agg = [results_agg[i][end] for i in 1:length(γ_array)]

end_results_hv = hcat(end_results_hv...)
end_results_av = hcat(end_results_av...)
end_results_agg = hcat(end_results_agg...)

edge_costs = travel_times.(end_results_agg, a, b)

costs = sum(end_results_agg .* edge_costs, dims=1)[:]
