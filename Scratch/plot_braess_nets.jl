using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      LinearAlgebra,
      SparseArrays

adj_mat = [0 1 1 0;
           0 0 1 1;
           0 0 0 1;
           0 0 0 0]

a = [0.5, 1.0, 0.1, 1.0, 0.5]
b = [1.0, 0.5, 0.1, 0.5, 1.0]

g_directed = DiGraph(adj_mat)

positions = [0.05 0.5;
             0.5 0.75;
             0.5 0.25;
             0.95 0.5]

rn = RoadNetwork(g_directed,
                 Dict(:a=>a, :b=>b), 
                 Dict(:pos=>positions))

n = nv(rn.g)
m = ne(rn.g)

# OD and penetration rate
# Note: strange arrays because its multi-od but we have only 1
ods = [(1,4)]
n_ods = length(ods)
d = 1.0
demands = [d]

#γ_array = [i/10 for i in 0:10]
γ_array = [0.1, 0.3, 0.5, 0.7, 0.9]

#= 
In this simple example with a sinlgle OD
we don't really need a tree and actually
using paths makes for easier examples.
=#

sorted_edges = collect(edges(rn.g))

# Names started making sense, then I had to add another one
good_path = [1, 4]
bad_path = [1, 3, 5]
double_path = [1, 4, 2, 5]


# This can all be shortened to one for loop...

# 3 cases of protected edges (store protected and complement)
# Case 1:
pe_1 = sorted_edges[good_path]
oi_1 = [] # other indices: Only AVs allowed
oe_1 = [] # other edges: Only AVs allowed
# Case 2:
pe_2 = sorted_edges[bad_path]
oi_2 = [] 
oe_2 = []
# Case 3
pe_3 = sorted_edges[double_path]
oi_3 = [] 
oe_3 = []

for (i, ed) in enumerate(sorted_edges)
    # Case 1
    if !in(ed, pe_1)
        push!(oi_1, i)
        push!(oe_1, ed)
    end
    # Case 2
    if !in(ed, pe_2)
        push!(oi_2, i)
        push!(oe_2, ed)
    end  
    # Case 3
    if !in(ed, pe_3)
        push!(oi_3, i)
        push!(oe_3, ed)
    end   
end

###
### Plotting of networks
###

#braess net with good path
ee = [(e.src, e.dst) for e in pe_1]
save_graph_tikz_edg(rn.g,
                    [ee],
                    positions,
                    "braess_p1.tex",
                    imp_nodes=[ods[1][1], ods[1][2]],
                    imp_labels=["O", "D"])
                    
#braess net with bad path
ee = [(e.src, e.dst) for e in pe_2]
save_graph_tikz_edg(rn.g,
                    [ee],
                    positions,
                    "braess_p2.tex",
                    imp_nodes=[ods[1][1], ods[1][2]],
                    imp_labels=["O", "D"])
                    
#braess net with double path
ee = [(e.src, e.dst) for e in pe_3]
save_graph_tikz_edg(rn.g,
                    [ee],
                    positions,
                    "braess_p3.tex",
                    imp_nodes=[ods[1][1], ods[1][2]],
                    imp_labels=["O", "D"])
