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
