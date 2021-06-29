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

###
### 
###

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

p_steps = [i for i in 0:length(r_array)]
γ_inds = [1, 3, 5]
for i in γ_inds

    hv_flows_bad = hcat(bad_results[1][i]...)
    hv_inner_bad = hv_flows_bad[3,:]
    hv_outer_bad = hv_flows_bad[2,:] + hv_flows_bad[4,:]
    
    av_flows_bad = hcat(bad_results[2][i]...)
    av_inner_bad = av_flows_bad[3,:]
    av_outer_bad = av_flows_bad[2,:] + av_flows_bad[4,:]
    
    plot(p_steps,
    [hv_inner_bad hv_outer_bad av_inner_bad av_outer_bad],
    labels = ["HVs inner" "HVs outer" "AVs inner" "AVs outer"])
    ylims!(0,d)
    
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Flow", yguidefontsize=16)
    title!("d=$(d), γ=$(γ_array[i])", titlefontsize=16)
    savefig("braess_av_hv_flows_d_$(d)_gamma_$(γ_array[i]).pdf")
end

###
### Plots of flows with control measures
###

### HVs forced on inner route

d_range = LinRange(0.001, 1.7, 150)
for i in γ_inds
    # HVs
    inner_flows_hv = zeros(length(d_range))
    outer_flows_hv = zeros(length(d_range))
    # AVs
    inner_flows_av = zeros(length(d_range))
    outer_flows_av = zeros(length(d_range))

    for (j,d) in enumerate(d_range)

        bad_results = me_excluded_assignment_penalty_pr(rn,
                                                         ods,
                                                         [d],
                                                         γ_array,
                                                         oi_2,
                                                         r_array)
        
        hv_flows_bad_converged = bad_results[1][i][end] 
        hv_inner_bad_converged = hv_flows_bad_converged[3]
        hv_outer_bad_converged = hv_flows_bad_converged[2] +
                                     hv_flows_bad_converged[4]
        
        av_flows_bad_converged = bad_results[2][i][end]
        av_inner_bad_converged = av_flows_bad_converged[3]
        av_outer_bad_converged = av_flows_bad_converged[2] + 
                                     av_flows_bad_converged[4]
        
        inner_flows_hv[j] = hv_inner_bad_converged
        outer_flows_hv[j] = hv_outer_bad_converged
        
        inner_flows_av[j] = av_inner_bad_converged
        outer_flows_av[j] = av_outer_bad_converged
    end

    agg_inner_flows = inner_flows_hv + inner_flows_av 
    agg_outer_flows = outer_flows_hv + outer_flows_av

    series = hcat(inner_flows_hv,
                  inner_flows_av,
                  agg_inner_flows,
                  outer_flows_hv,
                  outer_flows_av,
                  agg_outer_flows)
         
    labs = hcat("Inner HVs",
                "Inner AVs",
                "Total inner",
                "Outer HVs",
                "Outer AVs",
                "Total outer")

    plot(d_range,
         series,
         labels=labs)
    
    ylims!(0, d_range[end])
    
    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Flow", yguidefontsize=16)
    title!("HVs on inner route (γ=$(γ_array[i]))", titlefontsize=16)
    
    savefig("braess_flows_restricted_to_inner_gamma_$(γ_array[i]).pdf")
end

### HVs forced on outer routes

for i in γ_inds
    # HVs
    inner_flows_hv = zeros(length(d_range))
    outer_flows_hv = zeros(length(d_range))
    # AVs
    inner_flows_av = zeros(length(d_range))
    outer_flows_av = zeros(length(d_range))

    for (j,d) in enumerate(d_range)

        double_results = me_excluded_assignment_penalty_pr(rn,
                                                         ods,
                                                         [d],
                                                         γ_array,
                                                         oi_3,
                                                         r_array)
        
        hv_flows_double_converged = double_results[1][i][end] 
        hv_inner_double_converged = hv_flows_double_converged[3]
        hv_outer_double_converged = hv_flows_double_converged[2] +
                                     hv_flows_double_converged[4]
        
        av_flows_double_converged = double_results[2][i][end]
        av_inner_double_converged = av_flows_double_converged[3]
        av_outer_double_converged = av_flows_double_converged[2] + 
                                     av_flows_double_converged[4]
        
        inner_flows_hv[j] = hv_inner_double_converged
        outer_flows_hv[j] = hv_outer_double_converged
        
        inner_flows_av[j] = av_inner_double_converged
        outer_flows_av[j] = av_outer_double_converged
    end

    agg_inner_flows = inner_flows_hv + inner_flows_av 
    agg_outer_flows = outer_flows_hv + outer_flows_av

    series = hcat(inner_flows_hv,
                  inner_flows_av,
                  agg_inner_flows,
                  outer_flows_hv,
                  outer_flows_av,
                  agg_outer_flows)
         
    labs = hcat("Inner HVs",
                "Inner AVs",
                "Total inner",
                "Outer HVs",
                "Outer AVs",
                "Total outer")

    plot(d_range,
         series,
         labels=labs)
    
    ylims!(0, d_range[end])
    
    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Flow", yguidefontsize=16)
    title!("HVs on outer routes (γ=$(γ_array[i]))", titlefontsize=16)
    
    savefig("braess_flows_restricted_to_outers_gamma_$(γ_array[i]).pdf")
end
