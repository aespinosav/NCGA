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

    ###
    ###
    ###

    #=
      The following containers are 2D arrays,
      columns: γ
      rows: penalty steps
    =#

    # Running mean and std containers (Welford alorithm)
    r_mean_tot = zeros(penalty_steps + 1, length(γ_array))
    r_std_tot = zeros(penalty_steps + 1, length(γ_array))

    r_mean_normed_tot_costs = zeros(penalty_steps + 1, length(γ_array))
    r_std_normed_tot_costs = zeros(penalty_steps + 1, length(γ_array))

    r_mean_perv_hv = zeros(penalty_steps + 1, length(γ_array))
    r_std_perv_hv = zeros(penalty_steps + 1, length(γ_array))

    r_mean_perv_av = zeros(penalty_steps + 1, length(γ_array))
    r_std_perv_av = zeros(penalty_steps + 1, length(γ_array))

    r_mean_cost_diff = zeros(penalty_steps + 1, length(γ_array))
    r_std_cost_diff = zeros(penalty_steps + 1, length(γ_array))

    r_mean_percentage_links_hv = 0
    r_std_percentage_links_hv = 0

    for (i,file) in enumerate(net_files)

        # load net

        net = loadgraph(file, MGFormat()) # Undirected
        rn = skel2rn(net) # Directed

        # For simplicity of code below
        a = rn.edge_params[:a]
        b = rn.edge_params[:b]

        sorted_edges = collect(edges(rn.g))

        O = node_closest_to(net, [0,0]) # lower left-most node
        D = node_closest_to(net, [1,1]) # upper right-most node

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

        # Update mean of % of network dedicated to HVs
        hv_percent_edges = length(hv_permitted_edges) / ne(rn.g)
        r_mean_percentage_links_hv,  r_std_percentage_links_hv = update_sums_welford(i, r_mean_percentage_links_hv , r_std_percentage_links_hv, hv_percent_edges)



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
        =#

        for (j,γ) in enumerate(γ_array)

            # for each γⱼ (jth value of γ_array)
            flows_hv = hcat(results[1][j]...)
            flows_av = hcat(results[2][j]...)
            flows_agg = hcat(results[3][j]...)

            ####################################

            edge_costs = travel_times.(flows_agg[j], a, b)

            # Total costs (and normalisation)
            tot_costs = mapslices(x -> total_cost(x, a, b), flows_agg, dims=1)[:]
            normed_tot_costs = (tot_costs .- so_cost) ./ (ue_cost - so_cost)

            # Per-veh costs (and difference)
            perv_costs_hv = (sum((flows_hv .* edge_costs), dims=1) / d*(1-γ))[:]
            perv_costs_av = (sum((flows_av .* edge_costs), dims=1) / d*γ)[:]

            perv_cost_diff = perv_costs_hv - perv_costs_av

            # Welford algorithm for running std (mean included)
            r_mean_tot[:,j], r_std_tot[:,j] = update_sums_welford(i, r_mean_tot[:,j], r_std_tot[:,j], tot_costs)
            r_mean_normed_tot_costs[:,j], r_std_normed_tot_costs[:,j] = update_sums_welford(i,
                                                                                            r_mean_normed_tot_costs[:,j],
                                                                                            r_std_normed_tot_costs[:,j],
                                                                                            normed_tot_costs)

            r_mean_perv_hv[:,j], r_std_perv_hv[:,j] = update_sums_welford(i, r_mean_perv_hv[:,j], r_std_perv_hv[:,j], perv_costs_hv)
            r_mean_perv_av[:,j], r_std_perv_av[:,j] = update_sums_welford(i, r_mean_perv_av[:,j], r_std_perv_av[:,j], perv_costs_av)

            r_mean_cost_diff[:,j], r_std_cost_diff[:,j] = update_sums_welford(i, r_mean_cost_diff[:,j], r_std_cost_diff[:,j], perv_cost_diff)

        end


    end

    ### Welford
    std_tot_costs = sqrt.(calc_var_welford(N, r_std_tot))
    std_perv_costs_hv = sqrt.(calc_var_welford(N, r_std_perv_hv))
    std_perv_costs_av = sqrt.(calc_var_welford(N, r_std_perv_av))

    # % of usable links for HVs
    std_hv_link_percent = sqrt.(calc_var_welford(r_std_percentage_links_hv))

    # Normalised costs (relative to UE and SO)
    std_normed_tot_costs = sqrt.(calc_var_welford(N, r_std_normed_tot_costs))

    # Difference of per-veh costs
    std_perv_cost_diff = sqrt.(calc_var_welford(N, r_std_cost_diff))

    ###
    ### Write data files
    ###

    # Raw values
    writedlm("mean_total_costs.dat", r_mean_tot)
    writedlm("std_total_costs.dat", std_tot_costs)

    writedlm("mean_perv_costs_hv.dat", r_mean_perv_hv)
    writedlm("std_perv_costs_hv.dat", std_perv_costs_hv)

    writedlm("mean_perv_costs_av.dat", r_mean_perv_av)
    writedlm("std_perv_costs_av.dat", std_perv_costs_av)

    # Normalised
    writedlm("normed_mean_total_costs.dat", r_mean_normed_tot_costs)
    writedlm("std_normed_total_costs.dat", std_normed_tot_costs)

    # Per-vehicle cost difference
    writedlm("perv_mean_costs_diff.dat", r_mean_cost_diff)
    writedlm("perv_std_cost_diff.dat", std_perv_cost_diff)

    # Write parameter file
    save("batch_run_params.jld2", parsed_args)
end

main()
