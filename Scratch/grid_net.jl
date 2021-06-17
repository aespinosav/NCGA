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
      Ipopt,
      Plots

      
###
### Gurobi licence only once
###      

      
if !(@isdefined env)
    const env = Gurobi.Env()
end


###
### Incidence matrix correction (and other functions)
###


include("tools.jl")
include("penalty.jl")


###
### γ range and demand
###


γ_array = [0.1, 0.5, 0.9]
#d = 1.0
d = 1.5


###
### Define graph and network and MST
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


###
### Demand structure
###


ods = [(1,9)]
n_ods = length(ods)
demands = [d]


###
### Restricting edges for HVs
###


# Find minimum spanning tree

mst = kruskal_mst(g)
mst_rev = reverse.(mst)

mst_edges = [(e.src, e.dst) for e in mst]
mst_rev_edges = [(e.src, e.dst) for e in mst_rev]

# Identify MST edges
mst_indices = findfirst.(isequal.(mst), [sorted_edges])
mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])


###
### Protected indices and edges (MST)
###


# MST
pi_1 = vcat(mst_indices, mst_rev_indices)
pe_1 = vcat(mst, mst_rev)

# Other tree
pi_2 = sort([2, 8, 10, 18, 19, 21, 20, 14, 22, 24, 17, 23, 15, 7, 6, 4])
pe_2 = sorted_edges[pi_2]

# Edges not in MST
oi_1 = [] 
oe_1 = []
# Edges not in other tree
oi_2 = [] 
oe_2 = []
for (i, ed) in enumerate(sorted_edges)
    # MST
    if !in(ed, pe_1)
        push!(oi_1, i)
        push!(oe_1, ed)
    end
    # Other tree
    if !in(ed, pe_2)
        push!(oi_2, i)
        push!(oe_2, ed)
    end  
end


###
### Penalty STAP (Mixed assignment and excluded roads)
###

# Penalty parameter array
r_array = [10^i for i in 1:10]


mst_results = me_excluded_assignment_penalty_pr(rn,
                                                ods,
                                                demands,
                                                γ_array,
                                                oi_1,
                                                r_array)

tr2_results = me_excluded_assignment_penalty_pr(rn,
                                                ods,
                                                demands,
                                                γ_array,
                                                oi_2,
                                                r_array)
                                                

###
### Costs and Plots
###

for (k,γ) in enumerate(γ_array)
    
    # Total and per-Veh costs (MST)
    TC_array = [total_cost(mst_results[3][k][i], a, b)
                for i in 1:length(r_array)+1]
                
    HVC_array = [partial_cost(mst_results[1][k][i],
                              mst_results[3][k][i],
                              a,
                              b) for i in 1:length(r_array)+1] ./ (1-γ)*d
                 
    AVC_array = [partial_cost(mst_results[2][k][i],
                              mst_results[3][k][i],
                              a,
                              b) for i in 1:length(r_array)+1] ./ γ*d
    
    # Total and per-Veh costs (other tree)
    TC_array_t2 = [total_cost(tr2_results[3][k][i], a, b)
                   for i in 1:length(r_array)+1]
                   
    HVC_array_t2 = [partial_cost(tr2_results[1][k][i],
                                 tr2_results[3][k][i],
                                 a,
                                 b) for i in 1:length(r_array)+1] ./ (1-γ)*d
                    
    AVC_array_t2 = [partial_cost(tr2_results[2][k][i],
                                 tr2_results[3][k][i],
                                 a,
                                 b) for i in 1:length(r_array)+1] ./ γ*d
                                 
    ###
    ### Plotting
    ###
    
    pen_iters = 0:length(r_array)                             
    
    # Total cost plot
    
    plt = plot(pen_iters,
               [TC_array TC_array_t2],
               label=["Cost MST" "Cost tree 2"],
               markershape=[:circle :square],
               markeralpha = 0.6,
               legend=:bottomright,
               color_palette=palette(:tab10))
               
    ylims!(5,20)
    
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Cost", yguidefontsize=16)
    title!("Total costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt, "total_costs_2_trees_gamma_$(γ)_tst.pdf")
    
    # Per-vehicle cost plot
    
    plt2 = plot(pen_iters,
            [HVC_array AVC_array],
            label=["HV cost MST" "AV cost MST"],
            markershape=[:dtriangle :circle],
            markeralpha = 0.6,
            color_palette=palette(:tab10)[[1,10]])
            
    plot!(pen_iters,
          [HVC_array_t2 AVC_array_t2],
          label=["HV cost tree 2" "AV cost tree 2"],
          markershape=[:dtriangle :circle],
          markeralpha = 0.6,
          color_palette=palette(:tab10)[[2,4]])
          
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Class cost (per vehicle)", yguidefontsize=16)
    title!("Per-vehicle costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt2, "per_veh_costs_2_trees_gamma_$(γ)_tst.pdf")    
end
