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
### Function definitions
###

include("tools.jl")
include("penalty.jl")

# Command line args parsing (in tools.jl)
parsed_args = parse_commandline_drange()

###
### Gurobi global varaible for licence token
###

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

    ### Arguments from command line:

    # Optimiser choice
    optimiser = parsed_args["optim"]

    # Penalty parameters
    penalty_steps = parsed_args["penalty_steps"]
    r_array = [10^i for i in 1:penalty_steps]

    # Penetration rate range
    γ_init = parsed_args["gamma_init"]
    γ_step = parsed_args["gamma_step"]
    γ_stop = parsed_args["gamma_stop"]
    # Penetration rate
    γ_array = γ_init:γ_step:γ_stop

    # Demand range
    d_start = parsed_args["d_start"]
    d_step =  parsed_args["d_step"]
    d_stop =  parsed_args["d_stop"]
    # d range
    demand_range = d_start:d_step:d_stop

    # MH parameters
    μ = parsed_args["μ"]
    p_splice = parsed_args["p_splice"]

    # perhaps should be a function of the network size
    # Maybe the same number of paths should be added?
    num_metropolis_steps = parsed_args["mh_steps"]
    stagger = parsed_args["path_sample_stagger"]

    # Dirs and net files
    ens_dir = parsed_args["dir"]
    nn = parsed_args["num_nets"]

    ############################################

    # Directories

    files_in_dir = readdir(ens_dir)
    network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
    net_files = sort(network_files, lt=natural)

    # Use number of nets given, or all of them
    if nn > 0 && nn < length(net_files)
        net_files = net_files[1:nn]
    else
        parsed_args["num_nets"] = length(net_files)
    end
    N = length(net_files)

    #############################################

    # Shorthand for readability
    ld = length(demand_range)
    lg = length(γ_array)

    # Running mean and std containers (Welford alorithm)
    r_mean_tot = zeros(ld, lg)
    r_std_tot = zeros(ld, lg)

    r_mean_normed_tot_costs = zeros(ld, lg)
    r_std_normed_tot_costs = zeros(ld, lg)

    r_mean_perv_hv = zeros(ld, lg)
    r_std_perv_hv = zeros(ld, lg)

    r_mean_perv_av = zeros(ld, lg)
    r_std_perv_av = zeros(ld, lg)

    # UE and SO costs
    r_mean_ue = zeros(ld)
    r_mean_so = zeros(ld)
    r_std_ue = zeros(ld)
    r_std_so = zeros(ld)

    #r_mean_cost_diff = zeros(length(γ_array), length(demand_array))
    #r_std_cost_diff = zeros(length(γ_array), length(demand_array))

    #r_mean_percentage_links_hv = 0
    #r_std_percentage_links_hv = 0

    for (i,file) in enumerate(net_files)

        # Load net
        net = loadgraph(file, MGFormat()) # Undirected
        rn = skel2rn(net) # Directed

        a = rn.edge_params[:a]
        b = rn.edge_params[:b]

        sorted_edges = collect(edges(rn.g))

        O = node_closest_to(net, [0,0]) # lower left-most node
        D = node_closest_to(net, [1,1]) # upper right-most node

        # Identify edges

        av_excl_edges = make_exclusive_edge_set(net,
                                                sorted_edges,
                                                O,
                                                D,
                                                μ,
                                                p_splice,
                                                weight_func,
                                                num_metropolis_steps,
                                                stagger)

        av_excl_edges_inds = findfirst.(isequal.(av_excl_edges), [sorted_edges])



        for (k,d) in enumerate(demand_range)

            ###
            ### Assignment results
            ###

            ### UE and SO

            # Note: These functions use gurobi by default
            ue_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:ue)
            so_flows = multi_pair_stap_nc(rn, [(O,D)], d, regime=:so)

            ue_cost = total_cost(ue_flows, a, b)
            so_cost = total_cost(so_flows, a, b)

            update_sums_welford!(i, view(r_mean_ue, k), view(r_std_ue, k), ue_cost)
            update_sums_welford!(i, view(r_mean_so, k), viev(r_std_so, k), so_cost)

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

            for (j, γ) in enumerate(γ_array)

                ###
                ### Flows and edge costs
                ###

                # for each γⱼ (jth value of γ_array)
                flows_hv = results[1][j][end]
                flows_av = results[2][j][end]
                flows_agg = results[3][j][end]

                # Edge costs
                edge_costs = travel_times(flows_agg, a, b)


                ###
                ### Ensemble calculations
                ###

                # Total costs (and normalisation)
                tot_cost = edge_costs ⋅ flows_agg
                normed_tot_cost = (tot_cost - so_cost) / (ue_cost - so_cost)

                # Per-vehicle costs
                perv_costs_hv = (flows_hv ⋅ edge_costs) / (d*(1-γ))
                perv_costs_av = (flows_av ⋅ edge_costs) / (d*γ)


                ### Update containers

                # Total costs
                update_sums_welford!(i,
                                    view(r_mean_tot, k, j),
                                    view(r_std_tot, k, j),
                                    tot_cost)

                # Normalised costs
                update_sums_welford!(i,
                                    view(r_mean_normed_tot_costs, k, j),
                                    view(r_std_normed_tot_costs, k, j),
                                    normed_tot_cost)

                # Per-vehicle costs
                update_sums_welford!(i,
                                    view(r_mean_perv_hv, k, j),
                                    view(r_std_perv_hv, k, j),
                                    perv_costs_hv)

                update_sums_welford!(i,
                                    view(r_mean_perv_av, k, j),
                                    view(r_std_perv_av, k,j),
                                    perv_costs_av)




                # Total costs (and normalisation)
                #tot_costs = mapslices(x -> total_cost(x, a, b), flows_agg, dims=1)[:]
                #normed_tot_costs = (tot_costs .- so_cost) ./ (ue_cost - so_cost)

                # Per-veh costs (and difference)
                #perv_costs_hv = (sum((flows_hv .* edge_costs), dims=1) / (d*(1-γ)))[:]
                #perv_costs_av = (sum((flows_av .* edge_costs), dims=1) / (d*γ))[:]

                #perv_cost_diff = perv_costs_hv - perv_costs_av

                # Update containers
                #r_mean_tot[:,j], r_std_tot[:,j] = update_sums_welford(i, r_mean_tot[:,j], r_std_tot[:,j], tot_costs)

            end  # γ loop
        end  # d loop
    end  # net loop

    # Welford eval std
    r_std_tot= sqrt.(calc_var_welford(N, r_std_tot))
    r_std_normed_tot_costs = sqrt.(calc_var_welford(N, r_std_normed_tot_costs))
    r_std_perv_hv = sqrt.(calc_var_welford(N, r_std_perv_hv))
    r_std_perv_av = sqrt.(calc_var_welford(N, r_std_perv_av))

    r_std_ue = sqrt.(calc_var_welford(N, r_std_ue))
    r_std_so = sqrt.(calc_var_welford(N, r_std_so))

    ###
    ### Save files
    ###

    # Total cost
    writedlm("d_range_mean_tot_costs.dat", r_mean_tot)
    println("file saved: d_range_mean_tot_costs.dat")
    writedlm("d_range_std_tot_costs.dat", r_std_tot)
    println("file saved: d_range_std_tot_costs.dat")

    # Normalised cost
    writedlm("d_range_mean_normed_tot_costs.dat", r_mean_normed_tot_costs)
    println("file saved: d_range_mean_normed_tot_costs.dat")
    writedlm("d_range_std_normed_tot_costs.dat", r_std_normed_tot_costs)
    println("file saved: d_range_std_normed_tot_costs.dat")

    # Per-vehicle cost
    writedlm("d_range_mean_perv_costs_hv.dat", r_mean_perv_hv)
    println("file saved: d_range_mean_perv_costs_hv.dat")
    writedlm("d_range_std_perv_costs_hv.dat", r_std_perv_hv)
    println("file saved: d_range_std_perv_costs_hv.dat")

    writedlm("d_range_mean_perv_costs_av.dat", r_mean_perv_av)
    println("file saved: d_range_mean_perv_costs_av.dat")
    writedlm("d_range_std_perv_costs_av.dat", r_std_perv_av)
    println("file saved: d_range_std_perv_costs_av.dat")

    # UE and SO costs (col1:mean , col2:std)
    writedlm("d_range_mean_ue_costs.dat", hcat(r_mean_ue, r_std_ue))
    println("saved file: d_range_ue_costs.dat")
    writedlm("d_range_mean_so_costs.dat", hcat(r_mean_so, r_std_so))
    println("saved file: d_range_so_costs.dat")

    save("batch_run_params_d_range.jld2", parsed_args)
end

# Run main function
main()
