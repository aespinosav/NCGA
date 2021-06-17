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
### Cost functions
###

include("tools.jl")

###
### Braess network definition
###

adj_mat = [0 1 1 0;
           0 0 1 1;
           0 0 0 1;
           0 0 0 0]

a = [0.5, 1.0, 0.1, 1.0, 0.5]
b = [1.0, 0.5, 0.1, 0.5, 1.0]

positions = [0.05 0.5;
             0.5 0.75;
             0.5 0.25;
             0.95 0.5]

g_directed = DiGraph(adj_mat)
rn = RoadNetwork(g_directed,
                 Dict(:a=>a, :b=>b), 
                 Dict(:pos=>positions))

###                 
### OD and demands
###

ods = [(1,4)] # 1 OD pair
n_ods = length(ods)

d_range = LinRange(0.001, 1.7, 100)

ue_flows = [multi_pair_stap_nc(rn, ods, d, regime=:ue) for d in d_range]
so_flows = [multi_pair_stap_nc(rn, ods, d, regime=:so) for d in d_range]

###
### Route flows
###

# Route 1 ([1, 4])
r1_ue_flows = [ue_flows[i][4] for i in 1:length(ue_flows)]
r1_so_flows = [so_flows[i][4] for i in 1:length(so_flows)]
# Route 2 ([2, 5])
r2_ue_flows = [ue_flows[i][2] for i in 1:length(ue_flows)]
r2_so_flows = [so_flows[i][2] for i in 1:length(so_flows)]

outer_ue_flows = r1_ue_flows + r2_ue_flows
outer_so_flows = r1_so_flows + r2_so_flows

# Route 3 ([1, 3, 5])
r3_ue_flows = [ue_flows[i][3] for i in 1:length(ue_flows)]
r3_so_flows = [so_flows[i][3] for i in 1:length(so_flows)]


###
### Plot (Braess inner and outer route flows)
###

# Plots backend (GR is giving problems with vline)
pyplot()

plot(d_range,
     [outer_ue_flows outer_so_flows r3_ue_flows r3_so_flows],
     linewidth = 2,
     label = ["Outer UE" "Outer SO" "Inner UE" "Inner SO"],
     legend = :topleft)

# Outer route activation     

vline!([0.182, 0.364], linestyle=:dash, color=:grey, label=nothing)
vline!([0.804, 1.608], linestyle=:dash, color=:grey, label=nothing)

xlabel!("Demand", xguidefontsize=16)
ylabel!("Flow", yguidefontsize=16)

title!("UE and SO flows on Braess network", titlefontsize=16)

savefig("braess_ue_so_flows.pdf")
