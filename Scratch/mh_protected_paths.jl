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


###
### Load road network
###

# network generated with net_gen_medium.jl (n=8, α_hat=0.8, β=1.5)
mg = loadgraph("g_test_medium.mg", MGFormat())
rn = skel2rn(mg)

sorted_edges = collect(edges(rn.g))


###
### Parameters
###


# STAP parameters
ods = [(1, 64)]
d = 1.5

γ_array = [0.1, 0.5, 0.9]

n_ods = length(ods)
demands = [d]


# MH parameters
μ = 4
p_splice = 0.75

num_metropolis_steps = 5000

stagger = 5


# Penalty parameter array
r_array = [10^i for i in 1:7]


###
### Restricting edges for HVs
###

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
# stagger in steps still ahs to be properly selcted (this is a good "guess")
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

save_paths_tikz(rn.g,
                indep_paths,
                get_node_pos(mg),
                "tree_and_mhsample_paths_net.tex",
                imp_nodes=ods[1],
                imp_labels=["O", "D"],
                scale=10,
                standalone_doc=true)
                
edges_for_plotting = [(e.src, e.dst) for e in pe_both]
save_graph_tikz_edg(rn.g,
                    [edges_for_plotting],
                    get_node_pos(mg),
                    "tree_and_mhsample_paths_net_2.tex",
                    edge_labels=false,
                    imp_nodes=[1, 64],
                    imp_labels=["O", "D"])
                
writedlm("test_edge_indexed_paths.dat", indep_paths_ind)


###
### Optimise
###

results = me_excluded_assignment_penalty_pr(rn,
                                            ods,
                                            demands,
                                            γ_array,
                                            oi,
                                            r_array)
                                            

for (k,γ) in enumerate(γ_array)
    
    # Total and per-Veh costs (MST)
    TC_array = [total_cost(results[3][k][i],
                rn.edge_params[:a],
                rn.edge_params[:b])
                for i in 1:length(r_array)+1]
                
    HVC_array = [partial_cost(results[1][k][i],
                              results[3][k][i],
                              rn.edge_params[:a],
                              rn.edge_params[:b])
                 for i in 1:length(r_array)+1] ./ (1-γ)*d
                 
    AVC_array = [partial_cost(results[2][k][i],
                              results[3][k][i],
                              rn.edge_params[:a],
                              rn.edge_params[:b])
                 for i in 1:length(r_array)+1] ./ γ*d
                                 
    ###
    ### Plotting
    ###
    
    pen_iters = 0:length(r_array)                             
    
    # Total cost plot
    
    plt = plot(pen_iters,
               TC_array,
               label="Total cost",
               markershape=[:circle],
               markeralpha = 0.6,
               legend=:bottomright,
               color_palette=palette(:tab10))
    
    ylims!(80,140)
    
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Cost", yguidefontsize=16)
    title!("Total costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt, "total_costs_mhpaths_gamma_$(γ).pdf")
    
    # Per-vehicle cost plot
    
    plt2 = plot(pen_iters,
            [HVC_array AVC_array],
            label=["HV cost MST" "AV cost MST"],
            markershape=[:dtriangle :circle],
            markeralpha = 0.6,
            color_palette=palette(:tab10)[[1,10]])
          
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Class cost (per vehicle)", yguidefontsize=16)
    title!("Per-vehicle costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt2, "per_veh_costs_mhpaths_gamma_$(γ).pdf")    
end
