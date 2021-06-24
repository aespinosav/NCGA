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

"""
Find MST of `g` and return array of edges and indices

Finds MST and also gets reversed edges of it.
"""
function find_mst(g, sorted_edges)

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

#import LightGraphs.incidence_matrix
#function incidence_matrix(G)
#    I = vcat([src(e) for e in edges(G)], [dst(e) for e in edges(G)])
#    J = vcat(collect(1:ne(G)), collect(1:ne(G)))
#    V = vcat(fill(-1, ne(G)), fill(1, ne(G)))
#    return sparse(I, J, V)
#end

###
### Verification
###

#function forbidden_flows(x, forbidden_links)
#    x[forbidden_links]
#end

###
### Plotting for comparison
###

function flow_compare_plot(x, y)
    plt = scatter(x, label="HV")
    scatter!(y, label="AV")
    plt
end
