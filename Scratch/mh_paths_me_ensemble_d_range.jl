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
        "--d_start"
            help = "Demand range start"
            required = true
            arg_type = Float64
        "--d_step"
            help = "Demand range step size"
            required = true
            arg_type = Float64 
        "--d_stop"
            help = "Demand range stop"
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

    files_in_dir = readdir(ens_dir)
    network_files = filter(s -> split(s, '.')[end] == "mg", files_in_dir)
    net_files = sort(network_files, lt=natural)

    if nn > 0 && nn < length(net_files)
        net_files = net_files[1:nn]
    else
        parsed_args["num_nets"] = length(net_files)
    end
    N = length(net_files)

    #############################################

    # Running mean and std containers (Welford alorithm)
    r_mean_tot = zeros(length(γ_array), length(demand_array))
    r_std_tot = zeros(length(γ_array), length(demand_array))

    r_mean_normed_tot_costs = zeros(length(γ_array), length(demand_array))
    r_std_normed_tot_costs = zeros(length(γ_array), length(demand_array))

    r_mean_perv_hv = zeros(length(γ_array), length(demand_array))
    r_std_perv_hv = zeros(length(γ_array), length(demand_array))

    r_mean_perv_av = zeros(length(γ_array), length(demand_array))
    r_std_perv_av = zeros(length(γ_array), length(demand_array))

    r_mean_cost_diff = zeros(length(γ_array), length(demand_array))
    r_std_cost_diff = zeros(length(γ_array), length(demand_array))

    r_mean_percentage_links_hv = 0
    r_std_percentage_links_hv = 0

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

        av_excl_edges_inds = findfirst.(isequal.(av_excl_edges) , sorted_edges)


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
        av_excl_edges = setdiff(sorted_edges, hv_permitted_edges)
        av_excl_edges_inds = setdiff(1:length(sorted_edges), hv_permitter_edges_inds)


        for (k,d) in demand_range

            ###
            ### Assignment results
            ###

            ### UE and SO

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

            # for each γⱼ (jth value of γ_array)
            flows_hv = hcat(results[1][j]...)
            flows_av = hcat(results[2][j]...)
            flows_agg = hcat(results[3][j]...)

            # edge costs

            edge_costs = travel_times.(flows_agg[j], a, b)

            # Total costs (and normalisation)
            tot_costs = mapslices(x -> total_cost(x, a, b), flows_agg, dims=1)[:]
            normed_tot_costs = (tot_costs .- so_cost) ./ (ue_cost - so_cost)

            # Per-veh costs (and difference)
            perv_costs_hv = (sum((flows_hv .* edge_costs), dims=1) / (d*(1-γ)))[:]
            perv_costs_av = (sum((flows_av .* edge_costs), dims=1) / (d*γ))[:]

            perv_cost_diff = perv_costs_hv - perv_costs_av


        end
    end

end
