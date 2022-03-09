###
### Utility for setting up STAP
###

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
Function for finding indices of AV links in network that
correspond to AV exclusive links from genome.

Usage:
------

    a : Individual
    genome_map: Genome -> Network dictionary
"""
function av_links(a::T, genome_map) where {T<:AbstractIndividual}
    active_alleles = findall(x->x==1, a.genome)
    av_links = [genome_map[i] for i in active_alleles]
end


###
### Functions for finding costs (for affine cost functions))
###

function travel_times(x, a, b)
    a .+ (x .* b)
end

function marginal_travel_times(x, a, b)
    # lᵐ(x) = l(x) + x*l'(x)
    a .+ 2*(x .+ b)
end

function total_cost(x, a, b)
    tt = travel_times(x, a, b)
    tt ⋅ x
end

function partial_cost(x, tot_flow, a, b)
    tt = travel_times(tot_flow, a, b)
    tt ⋅ x
end
