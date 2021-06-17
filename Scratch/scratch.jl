using Plots

include("penalty_method.jl")
include("tools.jl")

a = rn.edge_params[:a]
b = rn.edge_params[:b]

###
### Convergernce of penalty method
###

# Change in flows

# Change in costs

TC_array = [total_cost(agg[i], a, b) for i in 1:length(r_array)+1]
HVC_array = [partial_cost(xx[i], agg[i], a, b) for i in 1:length(r_array)+1]./(1-γ)
AVC_array = [partial_cost(yy[i], agg[i], a, b) for i in 1:length(r_array)+1]./γ

###
### Plots
###

#pen_iters = 0:length(r_array)

# Total cost plot
#plt = scatter(pen_iters, TC_array, label="Total cost", legend=:bottomright)
#xlabel!("Penalty method iterations")
#ylabel!("Cost")
#title!("Penetration rate γ=$(γ)")

#savefig(plt, "penalty_example_totcosts_10_seq_$(γ)_split.pdf")

# Class costs
#plt2 = scatter(pen_iters,
#               [HVC_array AVC_array],
#               labels=["HV cost" "AV cost"])
#xlabel!("Penalty method iterations")
#ylabel!("Class cost (per vehicle)")
#title!("Penetration rate γ=$(γ)")

#savefig(plt2, "penalty_example_classcosts_10_seq_$(γ)_split.pdf")


#save_graph_tikz_edg(g_directed,
#                    [mst_edges, mst_rev_edges],
#                    positions,
#                    "test_net_mst_labels.tex",
#                    imp_nodes=[ods[1][1], ods[1][2]],
#                    imp_labels=["O", "D"])

###
### Plots comparing 2 trees
###

# Total cost plot
#pen_iters = 0:length(r_array)
#plt = plot(pen_iters, [TC_array TC_array_t2], label=["Cost MST" "Cost tree 2"], markershape=:circle, legend=:bottomright, color_palette=palette(:tab10))
#ylims!(0,40)
#xlabel!("Penalty method iterations")
#ylabel!("Cost")
#title!("Penetration rate γ=$(γ)")
#savefig(plt, "total_costs_2_trees_gamma_$(γ).pdf")


# Class costs
pen_iters = 0:length(r_array)
plt2 = plot(pen_iters,
            [HVC_array AVC_array],
            label=["HV cost MST" "AV cost MST"],
            markershape=[:dtriangle :circle],
            color_palette=palette(:tab10)[[1,10]])
plot!(pen_iters,
      [HVC_array_t2 AVC_array_t2],
      label=["HV cost tree 2" "AV cost tree 2"],
      markershape=[:dtriangle :circle],
      color_palette=palette(:tab10)[[2,4]])
xlabel!("Penalty method iterations")
ylabel!("Class cost (per vehicle)")
title!("Penetration rate γ=$(γ)")
#savefig(plt2, "per_veh_costs_2_trees_gamma_$(γ).pdf")
