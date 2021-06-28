using ArgParse,
      TrafficNetworks2,
      SkeletonCities2,
      MHPaths,
      LightGraphs,
      MetaGraphs,
      SimpleWeightedGraphs,
      GraphIO,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      DelimitedFiles,
      FileIO,
      NaturalSort

###
### Parse args (CL interface)
###

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
         "--optim"
            help = "Which optimiser: `gurobi` or `ipopt`?"
            arg_type = Symbol
            default = :gurobi
        "--penalty_steps"
            help = "Steps in the 'penalty method'"
            arg_type = Int
            default = 7
        "--gamma_init"
            help = "Penetration rate start value" 
            arg_type = Float64
            default= 0.1
        "--gamma_step"
            help = "Penetration rate start value" 
            arg_type = Float64
            default= 0.4
        "--gamma_stop"
            help = "Penetration rate start value" 
            arg_type = Float64
            default= 0.9
        "--dir"
            help = "Directory of files with networks in an ensemble"
            arg_type = String
            default = "."
        "--num_nets"
            help = "Number of networks to use"
            arg_type = Int
            default = -1
        "--mh_steps"
            help = "Transitions to attempt in MH algorithm"
            arg_type = Int
            default = 10000
        "--path_sample_stagger"
            help = "Spacing between 'independent' samples of MH paths"
            arg_type = Int
            default = 10
        "d"
            help = "Demand"
            required = true
            arg_type = Float64
        "μ"
            help = "Temperature parameter for MH path sampling"
            required = true
            arg_type = Float64
        "p_splice"
            help = "Splicing prob for MH step transition: 0 < p_splice < 1"
            required = true
            arg_type = Float64
    end

    return parse_args(s)
end


###
### Incidence matrix correction (and other functions)
###

include("tools.jl")
include("penalty.jl")

parsed_args = parse_commandline()

if parsed_args["optim"] == :gurobi
    using Gurobi
    # Use single licence token
    if !(@isdefined env)
        const env = Gurobi.Env()
    end
elseif parsed_args["optim"] == :ipopt
    using Ipopt
end

###
### Main
###

function main()

    # Optimiser choice
    optimiser = parsed_args["optim"]

    # Penalty parameters
    penalty_steps = parsed_args["penalty_steps"]
    r_array = [10^i for i in 1:penalty_steps]

    γ_init = parsed_args["gamma_init"]
    γ_step = parsed_args["gamma_step"]
    γ_stop = parsed_args["gamma_stop"]
    # Penetration rate
    γ_array = γ_init:γ_step:γ_stop

    # MH parameters
    μ = parsed_args["μ"]
    d = parsed_args["d"]
    p_splice = parsed_args["p_splice"]

    # perhaps should be a function of the network size
    # Maybe the same number of paths should be added?
    num_metropolis_steps = parsed_args["mh_steps"]
    stagger = parsed_args["path_sample_stagger"]

    # Dirs and net files
    ens_dir = parsed_args["dir"]
    nn = parsed_args["num_nets"]

    ############################################

    #for testing:
    #ens_dir = "/u/49/espinoa5/unix/Documents/Projects/NetworkEnsembles/BoundedAlfaBeta/10x10/alphahat_0.75_beta_1.5"
    #file = "g_test_medium.mg"

    files_in_dir = readdir(ens_dir)
    network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
    net_files = sort(network_files, lt=natural)

    if nn > 0 && nn < length(net_files)
        net_files = net_files[1:nn]
    else
        parsed_args["num_nets"] = length(net_files)
    end

    N = length(net_files)

    #=
      The following containers are 2D arrays,
      columns: γ
      rows: pen steps
    =#

    suma_tot = zeros(penalty_steps + 1, length(γ_array))
    suma_cuad_tot = zeros(penalty_steps + 1, length(γ_array))

    suma_perv_hv = zeros(penalty_steps + 1, length(γ_array))
    suma_perv_av = zeros(penalty_steps + 1, length(γ_array))

    for (i,file) in enumerate(net_files)

        # load net

        net = loadgraph(file, MGFormat()) # Undirected
        rn = skel2rn(net) # Directed

        # For simplicity of code below
        a = rn.edge_params[:a]
        b = rn.edge_params[:b]

        sorted_edges = collect(edges(rn.g))

        O = node_closest_to(net, [0,0])
        D = node_closest_to(net, [1,1])

        ############################################

        ###
        ### Identify edges (for HVs and exclusive for AVs)
        ###

        #MST
        mst_inds, mst_edges = find_mst(net, sorted_edges)

        # MH Paths
        mhp_nodes, mhp_edges, mhp_inds = find_mh_paths(net,
                                                       sorted_edges,
                                                       O,
                                                       D,
                                                       μ,
                                                       p_splice,
                                                       weight_func,
                                                       num_metropolis_steps,
                                                       stagger)
        mhp_all_edges_ind = vcat(mhp_inds...)
        mhp_all_edges = vcat(mhp_edges...)

        # HV permitted edges
        hv_permitted_edges_inds = unique(vcat(mst_inds, mhp_all_edges_ind))
        hv_permitted_edges = unique(vcat(mst_edges, mhp_all_edges))

        # AV exclusive edges
        av_excl_edges_inds = []
        av_excl_edges = []
        for (j, ed) in enumerate(sorted_edges)
            # MST
            if !in(ed, hv_permitted_edges)
                push!(av_excl_edges_inds, j)
                push!(av_excl_edges, ed)
            end
        end

        ############################################

        ###
        ### Assignment results
        ###

        ###
        ### UE and SO
        ###

        # Note: These functions use gurobi by default
        ue_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:ue)
        so_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:so)

        ue_cost = total_cost(ue_flows, a, b)
        so_cost = total_cost(so_flows, a, b)

        ###
        ### Mixed EQ on restricted network
        ###
        
        results = me_excluded_assignment_penalty_pr(rn,
                                                    [(O,D)],
                                                    [d],
                                                    γ_array,
                                                    av_excl_edges_inds,
                                                    r_array,
                                                    solver=optimiser)
        #=
        results[1] has HV flows
        results[2] has AV flows
        results[3] has aggregate flows
        
        Each results[i] array has as many elements as 
        there are γ_values ( i.e. length(γ_array) ).
        
        To make things easier to work with we transform to 
        an array of 2D arrays:
        =#

        flows_hv = [hcat(results[1][j]...) for j in 1:length(γ_array)]
        flows_av = [hcat(results[2][j]...) for j in 1:length(γ_array)]
        flows_agg = [hcat(results[3][j]...) for j in 1:length(γ_array)]

        # Reshaping might affect performance (maybe put all this in sum for updating?)
        edge_costs = [travel_times.(flows_agg[j], a, b) 
                      for j in 1:length(γ_array)]

        tot_costs = [mapslices(x -> total_cost(x, a, b), flows_agg[j], dims=1)
                     for j in 1:length(γ_array)]
        tot_costs = map(x->reshape(x, length(x), ), tot_costs)

        perv_costs_hv = [sum((flows_hv[j] .* edge_costs[j]), dims=1) / d*(1-γ_array[j])
                         for j in 1:length(γ_array)]
        perv_costs_hv = map(x->reshape(x, length(x), ), perv_costs_hv)


        perv_costs_av = [sum((flows_av[j] .* edge_costs[j]), dims=1) / d*γ_array[j]
                         for j in 1:length(γ_array)]
        perv_costs_av = map(x->reshape(x, length(x), ), perv_costs_av)

        # Update running sums
        for j in 1:length(γ_array)
            suma_tot[:,j] += tot_costs[j]
            suma_cuad_tot[:,j] += tot_costs[j].^2

            suma_perv_hv[:,j] += perv_costs_hv[j]
            suma_perv_av[:,j] += perv_costs_av[j]

            # for normalised values
        end
    end

    # Calculate means
    mean_tot_costs = suma_tot ./ N
    std_tot_costs = sqrt.((1/(N-1)) * (suma_cuad_tot .- N*mean_tot_costs.^2))

    mean_perv_costs_hv = suma_perv_hv ./ N
    mean_perv_costs_av = suma_perv_av ./ N

    # Normalised costs (relative to UE and SO)

    # Difference of per-veh costs

    ###
    ### Write data files
    ###

    # Raw values
    writedlm("mean_total_costs.dat", mean_tot_costs)
    writedlm("std_total_costs.dat", std_tot_costs)

    writedlm("mean_perv_costs_hv.dat", mean_perv_costs_hv)
    writedlm("mean_perv_costs_av.dat", mean_perv_costs_av)

    # Normalised

    # Write parameter file
    save("batch_run_params.jld2", parsed_args)
end

main()
