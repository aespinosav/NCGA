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
#μ = 3
#p_splice = 0.75
#num_metropolis_steps = 5000
#stagger = 5

nn = 10

# Assignment parameters
demand_array = 0.0001:0.0005:0.03
#demand_array = 0.00001:0.000005:0.001


# Dir with network files

#ens_dir = "/u/49/espinoa5/unix/Documents/Projects/NetworkEnsembles/BoundedAlfaBeta/10x10/alphahat_0.75_beta_1.5"
#ens_dir = "/u/49/espinoa5/unix/Documents/Projects/NetworkControl/Scratch/TestData"

###
### Get list of network files
###

files_in_dir = readdir(ens_dir)
network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
net_files = sort(network_files, lt=natural)

if length(net_files) >= nn
    net_files = net_files[1:nn]
end

##################################################################

PoA_array = zeros(length(demand_array), length(net_files))
peak_poa = zeros(length(net_files))

for (k, f) in enumerate(net_files)

    ###
    ### Load network and find paths
    ###

    filename = joinpath(ens_dir, f)

    #rn, ods, ex_edges = load_net_and_find_edges(filename,
    #                                            μ,
    #                                            p_splice,
    #                                            num_metropolis_steps,
    #                                            stagger)


    mg = loadgraph(filename, MGFormat())
    rn = skel2rn(mg)
    sorted_edges = collect(edges(rn.g))
    O = node_closest_to(mg, [0,0]) # lower left-most node
    D = node_closest_to(mg, [1,1]) # upper right-most node
    ods = [(O, D)]

    # Nicer notation
    a = rn.edge_params[:a]
    b = rn.edge_params[:b]

    ###
    ### Assignments (UE and SO)
    ###

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

    PoA_array[:,k] = poa
    peak_poa[k] = demand_array[argmax(poa)]
end

ens_name = split(ens_dir, "/")[end]

plot(demand_array, PoA_array, legend=false, grid=false)
ylabel!("PoA", yguidefontsize=16)
xlabel!("Demand", xguidefontsize=16)
lens!([0, 0.015], [1, 1.029], inset = (1, bbox(0.5, 0.0, 0.4, 0.4)))
title!(ens_name)
savefig(ens_name*"d_calibration_poa.pdf")
