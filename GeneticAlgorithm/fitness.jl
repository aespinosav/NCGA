# Functions for calculating fitness and interfacing with STAP code

#=
using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      GraphIO,
      MetaGraphs,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Ipopt
=#

###
### Fitnesses for `d` and `γ` ranges
###

"""
Calculates fitness of individual.

This fitness function does a very coarse integration of the total costs
in (d, γ)-space for the given demand and penetration rate ranges.

Very naive method, but we can see if it works.
"""
function fitness!(a::Individual,
             rn::RoadNetwork,
             genome_link_dict,
             ods,
             d_range::AbstractVector,
             γ_range::AbstractVector)

    n_ods = length(ods)

    ld = length(d_range)
    lg = length(γ_range)
    num_gridpoints = ld * lg

    cost_array = zeros(ld, lg)

    for (i,q) in enumerate(d_range)
        for (j,γ) in enumerate(γ_range)

            demands = (q/n_ods) * ones(n_ods)

            true_alleles = findall(x->x==1, a.genome)
            bound_links = [genome_link_dict[k] for k in true_alleles]

            x, y = stap_wrapper_fit(rn,
                                    ods,
                                    demands,
                                    γ,
                                    bound_links)

            agg = x + y

            tc = total_cost(agg, rn.edge_params[:a], rn.edge_params[:b])
            tc_av = tc / q

            cost_array[i,j] = tc_av
        end
    end
    fit = 1 / (sum(cost_array) / num_gridpoints)
    a.fitness = fit
end

function fitness!(pop::Array{Individual,1},
              rn::RoadNetwork,
              genome_link_dict,
              ods,
              d_range_or_float,
              γ_range_or_float)

    for ind in pop
        fitness!(ind, rn, genome_link_dict, ods, d_range_or_float, γ_range_or_float)
    end
end


###
### Fitnesses for single (d, γ) pair
###

"""
Fitness function for sigle (d, γ) pair.

Fitness is 1/average_costs
"""
function fitness!(a::Individual,
             rn::RoadNetwork,
             genome_link_dict::Dict{Int64,Int64},
             ods::Vector{Tuple{Int64, Int64}},
             d::Float64,
             γ::Float64)

    true_alleles = findall(x->x==1, a.genome)
    bound_links = [genome_link_dict[k] for k in true_alleles]

    n_ods = length(ods)
    od_demands = (d/n_ods) * ones(n_ods)

    x, y = stap_wrapper_fit(rn,
                            ods,
                            od_demands,
                            γ,
                            bound_links)

    agg = x + y

    tc = total_cost(agg, rn.edge_params[:a], rn.edge_params[:b])
    tc_av = tc / d

    fitness = 1/tc_av
    a.fitness = fitness
end


###
### Fitness defined in terms of dispersion of AV costs
###

"""
Fitness function for single (d, γ) value pair based on dispersion
of AV costs.

Single OD for now...

Naive function that calculates all simple paths (with threshold)
between O and D to get distribution of AV costs

Usage:
-----
    dispersion_fitness!(indi::Individual,
                        rn::RoadNetwork,
                        genome_link_dict::Dict{Int64,Int64},
                        ods::Vector{Tuple{Int64, Int64}},
                        d::Float64,
                        γ::Float64;
                        τ_av_flow::Float64=0.001,
                        p_cutoff::Int64=30,
                        τ_cost::Float64=0.15)

    + τ_av_flow: threshold for considering link unused
                (proportion of AV demand)

    + p_cutoff: max length of simple paths to find

    + τ_cost: proportion of marginal cost gap (max - min) to use as
              threshold for selecting min marginal cost paths
"""
function dispersion_fitness!(indi::Individual,
                             rn::RoadNetwork,
                             genome_link_dict::Dict{Int64,Int64},
                             ods::Vector{Tuple{Int64, Int64}},
                             d::Float64,
                             γ::Float64;
                             τ_av_flow::Float64=0.001,
                             p_cutoff::Int64=30,
                             τ_cost::Float64=0.15)

    ### Genome -> Network representation

    edgelist = collect(edges(rn.g))
    edgemap = Dict{Edge, Int}(e=>i for (i, e) in enumerate(edgelist))

    a = rn.edge_params[:a]
    b = rn.edge_params[:b]

    true_alleles = findall(x->x==1, indi.genome)
    bound_links = av_links(indi, genome_link_dict)

    ### STAP

    n_ods = length(ods)
    od_demands = (d/n_ods) * ones(n_ods)

    hv_flow, av_flow = stap_wrapper_fit(rn,
                                        ods,
                                        od_demands,
                                        γ,
                                        bound_links)
    agg_flow = hv_flow  + av_flow
    link_costs = travel_times(agg_flow, a, b)

    ### Sub-network selection

    threshold_av_flow = τ_av_flow * d * γ

    used_av_links = findall(x -> x>threshold_av_flow, av_flow)
    av_elist = edgelist[used_av_links]
    av_subgraph, av_vert_map = induced_subgraph(rn.g, av_elist)

    # Single OD (ods[1])
    av_O = findfirst(x->x==ods[1][1], av_vert_map)
    av_D = findfirst(x->x==ods[1][2], av_vert_map)

    ### Path selection

    # cutoff from argumetn p_cutoff
    path_iterator = all_simple_paths(av_subgraph,
                                     av_O, av_D, cutoff=p_cutoff)

    av_paths = Vector{Vector{Int64}}()
    for p in path_iterator
        Γ = edges_in_path(p, av_vert_map, edgemap)
        push!(av_paths, Γ)
    end

    ### Path costs

    # Path costs and marginal costs
    path_tts = Vector{Float64}(undef, length(av_paths))
    path_mar_tts = Vector{Float64}(undef, length(av_paths))
    for (i,p) in enumerate(av_paths)
        path_flow = agg_flow[p]
        path_a = a[p]
        path_b = b[p]

        path_cost = sum(travel_times(path_flow, path_a, path_b))
        mar_path_cost = sum(marginal_travel_times(path_flow, path_a, path_b))

        path_tts[i] = path_cost
        path_mar_tts[i] = mar_path_cost
    end

    # "Min marginal cost" paths

    max_mar_cost = maximum(path_mar_tts)
    min_mar_cost = minimum(path_mar_tts)
    # Min considered within 10% of spread
    threshold_mar_tt = min_mar_cost + τ_cost*(max_mar_cost - min_mar_cost)

    used_path_indices = findall(x->x<threshold_mar_tt, path_mar_tts)
    path_set = av_paths[used_path_indices]

    ### Determine path flows

    Λ = edge_route_incidence_matrix_ei(rn, path_set)

    # Find pseudoinverse (LSE) solution
    # A \ x is equivalent to pinv(A)*x (i.e. divide on the left)
    # Clip flows that become negative to 0

    ĥ = map(x -> x<0 ? 0.0  : x, Matrix(Λ) \ av_flow)

    # Consistent route flows with AV demand
    ĥ *= d*γ/sum(ĥ)

    # Get non-zero paths and flows
    non_z_inds = findall(x -> x>1e-12, ĥ)
    reduced_path_set = path_set[non_z_inds]
    # Redefine ĥ
    ĥ_r = ĥ[non_z_inds]
    Λ_r = Λ[:,non_z_inds]

    ### Path costs at estimated flows

    # Link flows under estimated path flows
    f̂ = Λ_r * ĥ_r
    new_agg_flow = f̂ + hv_flow
    ĉ = travel_times(new_agg_flow, a, b)

    # Path costs
    path_costs = Vector{Float64}(undef, length(reduced_path_set))
    for (i,p) in enumerate(reduced_path_set)
        p_cost = sum(travel_times(new_agg_flow[p], a[p], b[p]))
        path_costs[i] = p_cost
    end

    ### Dispersion measure and fitness

    # STD
    μ = mean(path_costs, StatsBase.weights(ĥ_r))
    σ = std(path_costs, StatsBase.weights(ĥ_r))

    # Fitenss inverse of coefficient of variation
    indi.fitness = μ/σ
end


"""
dispersion_fitness! applied to a population
"""
function dispersion_fitness!(pop::Array{Individual,1},
                             rn::RoadNetwork,
                             genome_link_dict,
                             ods,
                             d_range_or_float,
                             γ_range_or_float;
                             τ_av_flow::Float64=0.001,
                             p_cutoff::Int64=30,
                             τ_cost::Float64=0.15)
    for ind in pop
        dispersion_fitness!(ind,
                            rn,
                            genome_link_dict,
                            ods,
                            d_range_or_float,
                            γ_range_or_float,
                            τ_av_flow=τ_av_flow,
                            p_cutoff=p_cutoff,
                            τ_cost=τ_cost)
    end
end




### K-shortest paths version


"""
Fitness function for single (d, γ) value pair based on dispersion
of AV costs.

Single OD for now...

Naive function that calculates all simple paths (with threshold)
between O and D to get distribution of AV costs

Usage:
-----
    dispersion_fitness!(indi::Individual,
                        rn::RoadNetwork,
                        genome_link_dict::Dict{Int64,Int64},
                        ods::Vector{Tuple{Int64, Int64}},
                        d::Float64,
                        γ::Float64;
                        τ_av_flow::Float64=0.001,
                        p_cutoff::Int64=30,
                        τ_cost::Float64=0.15)

    + τ_av_flow: threshold for considering link unused
                (proportion of AV demand)

    + p_cutoff: max length of simple paths to find

    + τ_cost: proportion of marginal cost gap (max - min) to use as
              threshold for selecting min marginal cost paths
"""
function dispersion_fitness_ksp!(indi::Individual_KSP,
                             rn::RoadNetwork,
                             genome_link_dict::Dict{Int64,Int64},
                             ods::Vector{Tuple{Int64, Int64}},
                             d::Float64,
                             γ::Float64;
                             τ_av_flow::Float64=0.001,
                             k::Int64=20,
                             trimming::Bool=true)

    ### Genome -> Network representation

    edgelist = collect(edges(rn.g))
    edgemap = Dict{Edge, Int}(e=>i for (i, e) in enumerate(edgelist))

    a = rn.edge_params[:a]
    b = rn.edge_params[:b]

    bound_links = av_links(indi, genome_link_dict)

    ### STAP

    n_ods = length(ods)
    od_demands = (d/n_ods) * ones(n_ods)

    hv_flow, av_flow = stap_wrapper_fit(rn,
                                        ods,
                                        od_demands,
                                        γ,
                                        bound_links)
    agg_flow = hv_flow  + av_flow
    link_costs = travel_times(agg_flow, a, b)

    marginal_link_costs = marginal_travel_times(agg_flow, a, b)

    ### Sub-network selection

    threshold_av_flow = τ_av_flow * d * γ

    used_av_links = findall(x -> x>threshold_av_flow, av_flow)
    av_elist = edgelist[used_av_links]
    av_subgraph, av_vert_map = induced_subgraph(rn.g, av_elist)

    # Remove links that are selected as AV links but
    # carry no flow (forced mutation?), is this a good idea?
    if trimming == true
        remove_unused_av_links!(indi, genome_link_dict, av_flow; threshold=threshold_av_flow)
    end

    ### Path selection

    # Single OD (ods[1])
    av_O = findfirst(x->x==ods[1][1], av_vert_map)
    av_D = findfirst(x->x==ods[1][2], av_vert_map)

    av_paths = get_ksps(av_subgraph,
                        av_vert_map,
                        av_O,
                        av_D,
                        edgemap,
                        marginal_link_costs,
                        k)

    ### Path costs

    # Path costs and marginal costs
    path_tts = Vector{Float64}(undef, length(av_paths))
    path_mar_tts = Vector{Float64}(undef, length(av_paths))
    for (i,p) in enumerate(av_paths)
        path_flow = agg_flow[p]
        path_a = a[p]
        path_b = b[p]

        path_cost = sum(travel_times(path_flow, path_a, path_b))
        mar_path_cost = sum(marginal_travel_times(path_flow, path_a, path_b))

        path_tts[i] = path_cost
        path_mar_tts[i] = mar_path_cost
    end

    ### Determine path flows

    Λ = edge_route_incidence_matrix_ei(rn, av_paths)

    # Find pseudoinverse (LSE) solution and clip for >0 flows
    ĥ = map(x -> x<0 ? 0.0  : x, Matrix(Λ) \ av_flow)
    ĥ *= d*γ/sum(ĥ) # Consistent route flows with AV demand

    # Get non-zero paths and flows
    non_z_inds = findall(x -> x>1e-12, ĥ)
    reduced_path_set = av_paths[non_z_inds]

    ĥ_r = ĥ[non_z_inds]
    Λ_r = Λ[:,non_z_inds]

    ### Path costs at estimated flows

    # Link flows under estimated path flows
    f̂ = Λ_r * ĥ_r
    new_agg_flow = f̂ + hv_flow
    ĉ = travel_times(new_agg_flow, a, b)

    # Path costs
    path_costs = Vector{Float64}(undef, length(reduced_path_set))
    for (i,p) in enumerate(reduced_path_set)
        p_cost = sum(travel_times(new_agg_flow[p], a[p], b[p]))
        path_costs[i] = p_cost
    end

    ### Dispersion measure (coefficient of variation)
    μ = mean(path_costs, StatsBase.weights(ĥ_r))
    σ = std(path_costs, StatsBase.weights(ĥ_r))

    # Fitness: 1/CV

    indi.fitness = μ/σ

    # Add distribution details
    indi.mean_path_cost = μ
    indi.min_path_cost = minimum(path_costs)
    indi.max_path_cost = maximum(path_costs)
    indi.num_k_paths = length(path_costs)
end


function dispersion_fitness_ksp!(pop::Vector{Individual_KSP},
                             rn::RoadNetwork,
                             genome_link_dict::Dict{Int64,Int64},
                             ods::Vector{Tuple{Int64, Int64}},
                             d::Float64,
                             γ::Float64;
                             τ_av_flow::Float64=0.001,
                             k::Int64=20,
                             trimming::Bool=true)

    for indi in pop
        dispersion_fitness_ksp!(indi,
                                rn,
                                genome_link_dict,
                                ods,
                                d,
                                γ,
                                τ_av_flow=τ_av_flow,
                                k=k,
                                trimming=trimming)
    end
end

function dispersion_fitness!(pop::Vector{Individual_KSP},
                            rn::RoadNetwork,
                            genome_link_dict::Dict{Int64,Int64},
                            ods::Vector{Tuple{Int64, Int64}},
                            d::Float64,
                            γ::Float64;
                            τ_av_flow::Float64=0.001,
                            k::Int64=20,
                            trimming::Bool=true)

    dispersion_fitness_ksp!(pop,
                            rn,
                            genome_link_dict,
                            ods,
                            d,
                            γ,
                            τ_av_flow=τ_av_flow,
                            k=k,
                            trimming=trimming)
end


###
### Other functions
###


"""
Sorts a population based on fitness ranking.
Descending order
"""
function fitness_sort!(pop::Array{T,1}) where T <: AbstractIndividual
    sort!(pop, by= p -> p.fitness, rev=true)
end


"""
ods = [(O₁, D₁), ..., (Oₙ, Dₙ)]
demands: vector of demands for each OD pair
"""
function stap_wrapper_fit(rn, ods, demands, γ, bounded_indices; solver=:ipopt)

    # Structure
    n = nv(rn.g)
    m = ne(rn.g)
    A = incidence_matrix(rn.g)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])

    n_ods = length(ods)

    # Vectors for conservation constraints
    d_vects_hv = SparseVector{Float64,Int64}[]
    d_vects_av = SparseVector{Float64,Int64}[]
    for i in 1:length(ods)
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]

        push!(d_vects_hv, (1-γ).*d_vec)
        push!(d_vects_av, γ.*d_vec)
    end

    # Objective (Mixed multi-class Beckmann UE/SO)
    me_obj_func(x, y) = sum( B[i,i]*(x[i]+y[i])^2 +
                             a[i]*x[i] +
                             0.5*a[i]*y[i] for i in 1:length(x))

    # Choose solver
    if solver == :gurobi
        stap = Model(() -> Gurobi.Optimizer(env))
    elseif solver == :ipopt
        stap = Model(Ipopt.Optimizer)
    end

    # Supress output from solver
    set_silent(stap)

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)

    # OD specific link flows AV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    # Conservation constraints
    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])

    #Bound constraints (HV flow set to 0 on some links)
    for i in bounded_indices
        for k in 1:n_ods
            @constraint(stap, x[i,k]==0)
        end
    end

    # Unpenalised objective
    @objective(stap, Min, me_obj_func(x, y))

    # Solve
    optimize!(stap)
    return value.(x)[:], value.(y)[:]
end


###
### Functions for mapping AV paths back to edge indices
###


"""
Function to go from subgraph vertices back to
original graph vertices
"""
function back_to_verts(vert_list, vert_map)
    [vert_map[i] for i in vert_list]
end


"""
Get edge indices for AV paths
"""
function edges_in_path(av_sub_path, vert_map, edgemap)
    p_len = length(av_sub_path)

    oi = back_to_verts(av_sub_path, vert_map)
    edges_tups = [(oi[i], oi[i+1]) for i in 1:p_len-1]

    edges = [Edge(t[1], t[2]) for t in edges_tups]
    edge_indices = [edgemap[e] for e in edges]
end


"""
Returns Edge-Route incidence matrix (sparse) for a
given set of `paths` given as an array of arrays.
The paths shouls already be in terms of the edge indices.

Different from the function in TrafficNetwroks2
in the sense that it takes the edge indices as
ordered by edges(g).
"""
function edge_route_incidence_matrix_ei(rn, paths)
    m = ne(rn.g)
    r = length(paths)
    Λ = spzeros(Int, m, r)
    for (i,p) in enumerate(paths)
        Λ[p, i] .= 1
    end
    Λ
end


"""
Get k-shortest paths (using Yen's algorithm with
weights given by marginal_costs)
"""
function get_ksps(av_subgraph,
                  av_vert_map,
                  av_O,
                  av_D,
                  edgemap,
                  marginal_link_costs,
                  k)

    A = Graphs.adjacency_matrix(av_subgraph)
    I, J, V = findnz(A)
    num_entries = length(I)

    edge_marginal_tts = Vector{Float64}(undef, num_entries)
    @inbounds for i in 1:num_entries
        original_src = av_vert_map[I[i]]
        original_dst = av_vert_map[J[i]]
        ed_ind = edgemap[Edge(original_src, original_dst)]

        edge_marginal_tts[i] = marginal_link_costs[ed_ind]
    end

    A = sparse(I, J, edge_marginal_tts)
    yen_state = yen_k_shortest_paths(av_subgraph,
                                     av_O,
                                     av_D,
                                     A,
                                     k)

    av_paths = Vector{Vector{Int64}}()
    for p in yen_state.paths
        Γ = edges_in_path(p, av_vert_map, edgemap)
        push!(av_paths, Γ)
    end

    av_paths
end


"""
Removes AV paths if flow is below threshold,
the idea is to avoid surplus AV exclusive links
"""
function remove_unused_av_links!(indi, genome_map, av_flows; threshold=1e-6)
    for i in 1:length(indi.genome)
        if indi.genome[i] == 1
            if av_flows[genome_map[i]] < threshold
                indi.genome[i] = 0
            end
        end
    end
end
