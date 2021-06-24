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
### Def Braess network
###

adj_mat = [0 1 1 0;
           0 0 1 1;
           0 0 0 1;
           0 0 0 0]

a = [0.5, 1.0, 0.1, 1.0, 0.5]
b = [1.0, 0.5, 0.1, 0.5, 1.0]

g_directed = DiGraph(adj_mat)

positions = [0.05 0.5;
             0.5 0.75;
             0.5 0.25;
             0.95 0.5]

rn = RoadNetwork(g_directed,
                 Dict(:a=>a, :b=>b), 
                 Dict(:pos=>positions))

n = nv(rn.g)
m = ne(rn.g)

# OD and penetration rate
# Note: strange arrays because its multi-od but we have only 1
ods = [(1,4)]
n_ods = length(ods)
d = 1.0
demands = [d]

#γ_array = [i/10 for i in 0:10]
γ_array = [0.1, 0.3, 0.5, 0.7, 0.9]

#= 
In this simple example with a sinlgle OD
we don't really need a tree and actually
using paths makes for easier examples.
=#

sorted_edges = collect(edges(rn.g))

# Names started making sense, then I had to add another one
good_path = [1, 4]
bad_path = [1, 3, 5]
double_path = [1, 4, 2, 5]


# This can all be shortened to one for loop...

# 3 cases of protected edges (store protected and complement)
# Case 1:
pe_1 = sorted_edges[good_path]
oi_1 = [] # other indices: Only AVs allowed
oe_1 = [] # other edges: Only AVs allowed
# Case 2:
pe_2 = sorted_edges[bad_path]
oi_2 = [] 
oe_2 = []
# Case 3
pe_3 = sorted_edges[double_path]
oi_3 = [] 
oe_3 = []

for (i, ed) in enumerate(sorted_edges)
    # Case 1
    if !in(ed, pe_1)
        push!(oi_1, i)
        push!(oe_1, ed)
    end
    # Case 2
    if !in(ed, pe_2)
        push!(oi_2, i)
        push!(oe_2, ed)
    end  
    # Case 3
    if !in(ed, pe_3)
        push!(oi_3, i)
        push!(oe_3, ed)
    end   
end

# Penalty term array (starting from second iteration)
r_array = [10^i for i in 0:6]

# Actual assignment and optimisation
good_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_1,
                                                 r_array)

bad_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_2,
                                                 r_array)
                                                 
double_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_3,
                                                 r_array)
                                                 
                                                 
#=
The variables good_results and bad_results are arrays of arrays, 
their indices correspond to:
    
    1 - HV flows
    2 - AV flows
    3 - Aggregate flows
    
Each of these has n arrys as elements (one for each value of γ)

Each of these then has one array for each penalty iteration
=#

for γ_ind in [1, 3, 5] 

    # Case 1
    tot_cost = [total_cost(good_results[3][γ_ind][i], a, b) 
                for i in 1:length(r_array)]
                
    hv_cost_per_v = [partial_cost(good_results[1][γ_ind][i],
                                  good_results[3][γ_ind][i],
                                  a,
                                  b) 
                     for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v = [partial_cost(good_results[2][γ_ind][i],
                                  good_results[3][γ_ind][i],
                                  a,
                                  b)
                     for i in 1:length(r_array)] ./ γ_array[γ_ind]*d

    
    # Case 2
    tot_cost2 = [total_cost(bad_results[3][γ_ind][i], a, b) 
                 for i in 1:length(r_array)]
                 
                 
    hv_cost_per_v2 = [partial_cost(bad_results[1][γ_ind][i],
                                   bad_results[3][γ_ind][i],
                                   a,
                                   b) 
                      for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v2 = [partial_cost(bad_results[2][γ_ind][i],
                                   bad_results[3][γ_ind][i],
                                   a,
                                   b)
                      for i in 1:length(r_array)] ./ γ_array[γ_ind]*d
                     
    
    # Case 3
    tot_cost3 = [total_cost(double_results[3][γ_ind][i], a, b) 
                 for i in 1:length(r_array)]
    
    
    hv_cost_per_v3 = [partial_cost(double_results[1][γ_ind][i],
                                   double_results[3][γ_ind][i],
                                   a,
                                   b) 
                      for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v3 = [partial_cost(double_results[2][γ_ind][i],
                                   double_results[3][γ_ind][i],
                                   a,
                                   b)
                      for i in 1:length(r_array)] ./ γ_array[γ_ind]*d
    
    ###
    ### Plotting
    ###
    
    # Penalty iterations
    pen_iters = 0:length(r_array)-1
    
    # Total Costs (all 3 cases)
    # Legend labels
    l1 = "Cost (upper protected)"
    l2 = "Cost (interior protected)"
    l3 = "Cost (both outer routes protected)"
    
    plt = plot(pen_iters,
               [tot_cost tot_cost2 tot_cost3],
               label=[l1 l2 l3],
               linewidth = 2,
               markershape=[:utriangle :circle :square],
               markeralpha = 0.6,
               legend=:topleft,
               foreground_color_legend = nothing,
               background_color_legend = nothing,
               color_palette=palette(:tab10))
    
    #ylims!(0.7,0.75)
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Cost", yguidefontsize=16)
    title!("Total costs (d=$d, γ=$(γ_array[γ_ind]))", titlefontsize=16)
    
    savefig(plt, "total_costs_braess_gamma_d_$(d)_$(γ_array[γ_ind]).pdf")
    
   # Costs per vehicle (cases 2 and 3)   
   plt2 = plot(pen_iters,
          [hv_cost_per_v2 av_cost_per_v2],
          label=["HVs inner" "AVs inner"],
          linewidth = 2,
          markershape=[:dtriangle :hexagon],
          markeralpha = 0.6,
          legend=:bottomright,
          color_palette=palette(:tab10)[[1,10]])
         
    plot!(pen_iters,
         [hv_cost_per_v3 av_cost_per_v3],
         label=["HVs both outer" "AVs both outer"],
         linewidth = 2,
         markershape=[:utriangle :circle],
         markeralpha = 0.6,
         legend=:bottomright,
         foreground_color_legend = nothing,
         background_color_legend = nothing,
         color_palette=palette(:tab10)[[2,4]])
         
    ylabel!("Cost per vehicle", yguidefontsize=16)
    xlabel!("Penalty method iterations", xguidefontsize=16)
    title!("Per-vehicle costs (d=$d, γ=$(γ_array[γ_ind]))", titlefontsize=16)

    savefig(plt2, "perv_costs_braess_gamma_d_$(d)_$(γ_array[γ_ind]).pdf")
end
