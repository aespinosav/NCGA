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
      NaturalSort,
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
#filename = "g_test_medium.mg"

###
### Parameters
###

# MH parameters
μ = 4
p_splice = 0.75

num_metropolis_steps = 5000
stagger = 5

# Assignment parameters
demand_array = 0.0001:0.001:0.3

# Dir with network files

#ens_dir = "/u/49/espinoa5/unix/Documents/Projects/NetworkEnsembles/BoundedAlfaBeta/10x10/alphahat_0.75_beta_1.5"
ens_dir = "/u/49/espinoa5/unix/Documents/Projects/NetworkControl/Scratch/TestData"

###
### Get list of network files
###

files_in_dir = readdir(ens_dir)
network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
net_files = sort(network_files, lt=natural)


##################################################################

PoA_array = zeros(length(demand_array), length(net_files))
peak_poa = zeros(length(net_files))

for (k, f) in enumerate(net_files)

    ###
    ### Load network and find paths
    ###
    
    filename = joinpath(ens_dir, f)
    
    rn, ods, ex_edges = load_net_and_find_edges(filename,
                                                μ,
                                                p_splice,
                                                num_metropolis_steps,
                                                stagger)
    # Nicer notation
    a = rn.edge_params[:a]
    b = rn.edge_params[:b]

    ###
    ### Assignments (UE and SO)
    ###   

    UE_costs = zeros(length(demand_array))
    SO_costs = zeros(length(demand_array))

    for (k,d) in enumerate(demand_array)

        ue_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:ue)
        so_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:so)

        ue_cost = total_cost(ue_flows, a, b)
        so_cost = total_cost(so_flows, a, b)

        UE_costs[k] = ue_cost
        SO_costs[k] = so_cost
    end

    poa = UE_costs ./ SO_costs

    PoA_array[:,k] = poa
    peak_poa[k] = demand_array[argmax(poa)]
end

