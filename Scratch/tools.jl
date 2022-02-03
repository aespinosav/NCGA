###
### Parse args (CL interface)
###

using ArgParse,
      JuMP,
      LinearAlgebra

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


function parse_commandline_drange()
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
        "--num_path_samples"
            help = "If using different path sets, choose how many"
            arg_type = Int
            default = 5
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
### Import `multi_pair_stap_nc` and redefine to give solver options 
###

import TrafficNetworks2.multi_pair_stap_nc

"""
Generate and solve multi pair stap by making copies of
the network (general expression of stap)
"""
function multi_pair_stap_nc(rn, ods, demands; regime=:ue, solver=:gurobi)
    C = (regime == :ue ? 0.5 : 1.0)
    n = nv(rn.g)
    m = ne(rn.g)
    # Number of OD pairs
    n_ods = length(ods)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])

    d_vects = SparseVector{Float64,Int64}[]
    for i in 1:length(ods)
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]
        push!(d_vects, d_vec)
    end

    A = incidence_matrix(rn.g)

    if solver == :gurobi
        if @isdefined env
            stap = Model(() -> Gurobi.Optimizer(env))
        else
            stap = Model(Gurobi.Optimizer)
        end
    elseif solver == :ipopt
        stap = Model(Ipopt.Optimizer)
        #set_silent(stap)
    end

    #stap = Model(Gurobi.Optimizer)
    set_silent(stap)

    # OD specific link flows
    @variable(stap, x[1:m,1:n_ods] >= 0)
    # Aggregate link flows (for expressing objective)
    @variable(stap, link_flow[1:m] >= 0)

    @constraint(stap, [i in 1:n_ods], A*x[:,i] .== d_vects[i])
    # Link flows have to add up
    @constraint(stap,
                inter_var_con[i in 1:m], # name of constraint
                link_flow[i] == sum(x[i,j] for j in 1:n_ods))

    # Objective function
    #For some reason this stopped working...
    #@objective(stap, Min, dot(a,link_flow) + (C*link_flow'*B*link_flow))
    @objective(stap, Min, sum(a.*link_flow) + C*link_flow'*B*link_flow )
    optimize!(stap)
    value.(x)
end

###
### Update mean / variance (Welford)
###

function update_sums_welford(r_N, r_mean, r_M2, new_x)
    delta = new_x - r_mean
    r_mean += delta / r_N
    delta_2 = new_x - r_mean
    r_M2 += delta .* delta_2

    return r_mean, r_M2
end

#function update_sums_welford!(r_n::Real, r_mean::Real, r_m2::Real, new_x::Real)
#    delta = new_x - r_mean
#    r_mean += delta / r_n
#    delta_2 = new_x - r_mean
#    r_m2 += delta * delta_2
#end

function update_sums_welford!(r_n::Real, r_mean::AbstractArray, r_m2::AbstractArray, new_x)
    delta = new_x .- r_mean
    r_mean .+= delta / r_n
    delta_2 = new_x .- r_mean
    r_m2 .+= delta .* delta_2
end


# This function seems almost pointless...
"""
Running sample variance
"""
function calc_var_welford(N, M2)
    if N < 2
        return NaN
    else
        return M2/(N-1)
    end
end


###
### Cost functions
###

function travel_times(x, a, b)
    a .+ x.*b
end

function total_cost(x, a, b)
    tt = travel_times(x, a, b)
    tt ⋅ x
end

function partial_cost(x, tot_flow, a, b)
    tt = travel_times(tot_flow, a, b)
    tt ⋅ x
end

function my_incidence_matrix(G)
    I = vcat([src(e) for e in edges(G)], [dst(e) for e in edges(G)])
    J = vcat(collect(1:ne(G)), collect(1:ne(G)))
    V = vcat(fill(-1, ne(G)), fill(1, ne(G)))
    return sparse(I, J, V)
end

###
### Importing graph
###

"""
Loads graph from file (saved in MetaGraph format)
returns a road network

    `load_rn2(filename)`


"""
function load_rn2(filename)

    g = loadgraph(filename, MGFormat())
    positions = get_node_pos(g)

    # Make directed
    g_directed = DiGraph(g.graph)

    # Make road network
    a = edge_lengths(g_directed, positions)
    b = resource_allocation(g_directed, a)

    edge_params = Dict(:a => a, :b => b)
    node_params = Dict(:pos => positions)

    rn = RoadNetwork(g_directed, edge_params, node_params)
end


function weighted_adjacency(g)
    weighted_adj = Float64.(copy(adjacency_matrix(g)))
    for edge in edges(g)
        w = get_prop(g, edge, :length)
        weighted_adj[edge.src, edge.dst] = w
        weighted_adj[edge.dst, edge.src] = w
    end
    return weighted_adj
end

"""
Find node in network `mg` that is closes to the given `point`
"""
function node_closest_to(mg, point)
    positions = get_node_pos(mg)

    dists = zeros(nv(mg))
    min_dist = Inf
    min_node = -1
    for i in 1:nv(mg)
        dists[i] = norm(positions[i,:] - point)
        if dists[i] <= min_dist
           min_dist = dists[i]
           min_node = i
        end
    end
    return min_node
end

function node_closest_to(rn::RoadNetwork, point)
    positions = rn.node_params[:pos]
    num_nodes = nv(rn.g)

    dists = zeros(num_nodes)
    min_dist = Inf
    min_node = -1
    for i in 1:num_nodes
        dists[i] = norm(positions[i,:] - point)
        if dists[i] <= min_dist
           min_dist = dists[i]
           min_node = i
        end
    end
    return min_node
end

"""
Find MST of `mg` and return array of edges and indices

Finds MST and also gets reversed edges of it.
"""
function find_mst(mg, sorted_edges)

    #sorted_edges = collect(edges(mg))

    # Minimum Spanning Tree
    mst = kruskal_mst(mg)
    mst_rev = reverse.(mst)

    mst_edges = [(ed.src, ed.dst) for ed in mst]
    mst_rev_edges = [(ed.src, ed.dst) for ed in mst_rev]

    # Identify MST edges
    mst_indices = findfirst.(isequal.(mst), [sorted_edges])
    mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])

    # MST
    index_array = vcat(mst_indices, mst_rev_indices)
    edge_array = vcat(mst, mst_rev)

    return index_array, edge_array
end


"""
Sample paths with Metropolis-Hastings algorithm and return
arrays of paths expressed in terms of nodes, edges, and edge indices
"""
function find_mh_paths(mg,
                       sorted_edges,
                       O,
                       D,
                       μ,
                       p_splice,
                       weight_func,
                       num_metropolis_steps,
                       stagger)


    # Metropolis-Hastings (Flötteröd + Bierlaire)

    g_for_sampling = SimpleWeightedDiGraph(weighted_adjacency(mg))

    mh = MHInstance(g_for_sampling, O, D, μ, p_splice)
    mh_evolve!(mh, num_metropolis_steps, weight_func)

    paths = unique([p.Γ for p in mh.history])
    # stagger in steps still ahs to be properly selcted (this is a good "guess")
    indep_paths_nodes = paths[1:stagger:end]

    # paths are returned as lists of consecutive nodes

    indep_paths_edges = Array{Edge,1}[]
    for p in indep_paths_nodes
        ed_path = Array{Edge,1}(undef, length(p)-1)
        for (i, n) in enumerate(p[1:end-1])
            ed_path[i] = Edge(n, p[i+1])
        end
        push!(indep_paths_edges, ed_path)
    end

    indep_paths_inds = [findfirst.(isequal.(p), [sorted_edges])
                       for p in indep_paths_edges]

    return indep_paths_nodes, indep_paths_edges, indep_paths_inds
end


function make_exclusive_edge_set(net,
                                 sorted_edges,
                                 O,
                                 D,
                                 μ,
                                 p_splice,
                                 weight_func,
                                 num_metropolis_steps,
                                 stagger)

    # MST
    mst_inds, mst_edges = find_mst(net, sorted_edges)

    # Sample paths (Metropolis-Hastings)
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

    # Union of MST and MH paths
    hv_permitted_edges_inds = unique(vcat(mst_inds, mhp_all_edges_ind))
    hv_permitted_edges = unique(vcat(mst_edges, mhp_all_edges))

    # Edges for exclusive AV use
    av_excl_edges = setdiff(sorted_edges, hv_permitted_edges)

    return av_excl_edges
end

function make_exclusive_edge_set_ind(net,
                                    sorted_edges,
                                    O,
                                    D,
                                    μ,
                                    p_splice,
                                    weight_func,
                                    num_metropolis_steps,
                                    stagger)

    # MST
    mst_inds, mst_edges = find_mst(net, sorted_edges)

    # Sample paths (Metropolis-Hastings)
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

    # Union of MST and MH paths
    hv_permitted_edges_inds = unique(vcat(mst_inds, mhp_all_edges_ind))
    hv_permitted_edges = unique(vcat(mst_edges, mhp_all_edges))

    # Edges for exclusive AV use
    av_excl_edges = setdiff(sorted_edges, hv_permitted_edges)
    av_excl_edges_inds = findfirst.(isequal.(av_excl_edges), [sorted_edges])

    return av_excl_edges_inds
end

function excl_edge_set_sample(net,
                              sorted_edges,
                              O,
                              D,
                              μ,
                              p_splice,
                              weight_func,
                              num_metropolis_steps,
                              stagger,
                              samples)

    sampled_sets = Array{Int,1}[]
    for i in 1:samples
        s = make_exclusive_edge_set_ind(net,
                                        sorted_edges,
                                        O,
                                        D,
                                        μ,
                                        p_splice,
                                        weight_func,
                                        num_metropolis_steps,
                                        stagger)
        push!(sampled_sets, s)
    end

    return sampled_sets
end

###
### Bundle functions for easy importing of nets for this
###

"""
Bespoke function for experimenting on ensembles with this mixed equilibrium
and excluded edges for HVs.

Uses MH sampler

Returns:
    + rn :: RoadNetwork
    + ods :: Array of OD tuples
    + av_excl_edges_inds :: Array of indices of edges on which only AVs allowed
"""
function load_net_and_find_edges(filename,
                                 μ,
                                 p_splice,
                                 num_metropolis_steps,
                                 stagger;
                                 o_point=[0,0],
                                 d_point=[1,1])

    net = mg = loadgraph(filename, MGFormat())
    rn = skel2rn(mg)

    sorted_edges = collect(edges(rn.g))

    O = node_closest_to(net, o_point) # lower left-most node
    D = node_closest_to(net, d_point) # upper right-most node
    ods = [(O, D)]

    # Minimum Spanning Tree
    mst = kruskal_mst(mg)
    mst_rev = reverse.(mst)
    mst_edges = [(e.src, e.dst) for e in mst]
    mst_rev_edges = [(e.src, e.dst) for e in mst_rev]
    # Identify MST edges
    mst_indices = findfirst.(isequal.(mst), [sorted_edges])
    mst_rev_indices = findfirst.(isequal.(mst_rev), [sorted_edges])
    # MST
    pi_1 = vcat(mst_indices, mst_rev_indices)
    pe_1 = vcat(mst, mst_rev)

    # MH Paths
    mhp_nodes, mhp_edges, mhp_inds = find_mh_paths(net,
                                                   sorted_edges,
                                                   O,
                                                   D,
                                                   μ,
                                                   p_splice,
                                                   weight_func, # This is imported from module
                                                   num_metropolis_steps,
                                                   stagger)
    mhp_all_edges_ind = vcat(mhp_inds...)
    mhp_all_edges = vcat(mhp_edges...)
    # HV permitted edges
    hv_permitted_edges_inds = unique(vcat(mst_indices, mhp_all_edges_ind))
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

    return rn, ods, av_excl_edges_inds
end

###
### Plotting for comparison
###

function flow_compare_plot(x, y)
    plt = scatter(x, label="HV")
    scatter!(y, label="AV")
    plt
end
