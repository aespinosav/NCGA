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
### Test network
###

###
### Braess network
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
d = 0.4
demands = [d]

γ = 0.0

# Path protection
sorted_edges = collect(edges(rn.g))

double_path = [1, 4, 2, 5]
pe_3 = sorted_edges[double_path]
oi_3 = [] 
oe_3 = []
for (i, ed) in enumerate(sorted_edges)
    if !in(ed, pe_3)
        push!(oi_3, i)
        push!(oe_3, ed)
    end   
end


# Single penalty method iteration
r_array = [10^6]

double_results = me_excluded_assignment_penalty(rn,
                                                ods,
                                                demands, 
                                                γ,
                                                oi_3,
                                                r_array)
eq_flows = double_results[1][1]                                                
ue_flows = multi_pair_stap_nc(rn, ods, demands; regime=:ue)
