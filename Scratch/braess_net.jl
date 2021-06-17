using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO,
      StatsBase,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi,
      Plots

###
### Gurobi licence only once
###      
      
if !(@isdefined env)
    const env = Gurobi.Env()
end
      
###
### Incidence matrix correction (and other functions)
###

include("tools.jl")


###
### Function for assignment
###

function me_excluded_assignment_penalty(rn, ods, demands, γ, bounded_indices, r_array) 
    
    # Cost structure (standard)
    A = my_incidence_matrix(rn.g)
    a = rn.edge_params[:a]
    B = Diagonal(rn.edge_params[:b])
    
    # Vectors for conservation constraints (including OD terms)
    n_ods = length(ods)
    
    d_vects_hv = SparseVector{Float64,Int64}[]
    d_vects_av = SparseVector{Float64,Int64}[]
    for i in 1:n_ods
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]
        
        push!(d_vects_hv, (1-γ).*d_vec)
        push!(d_vects_av, γ.*d_vec)
    end
    
    # ME objective function (unpenalised)       
    me_obj_func(x,y) = ((x+y)'*B*(x+y) + a'*x + 0.5*a'*y)[1]
                             
    # Penalty function per link
    penalty_function(x, r, m) = r*(x^m)
    
    # penalty term for objective
    function pen_term(x, param, inds; power=2)
       ps = penalty_function.(x[inds], 0.5*a[inds], power)
       param*sum(ps)
    end
    
    # Penalised objective
    p_obj(x, y, p_param, b_indices) = me_obj_func(x,y) + 
                                      pen_term(x, p_param, b_indices)


    ### JuMP model construction
    
    #stap = Model(Ipopt.Optimizer)
    #stap = Model(Gurobi.Optimizer)
    stap = Model(() -> Gurobi.Optimizer(env))

    # OD specific link flows HV         
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, non_neg_x, x .>= 0)

    # OD specific link flows AV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, non_neg_y, y .>= 0)

    #hv_flow = sum(x, dims=2)
    #av_flow = sum(y, dims=2)

    # Conservation constraints
    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])


    # Unpenalised objective
    @objective(stap, Min, me_obj_func(x, y))
    
    ### 
    ### Solve unpenalised
    ###
    
    optimize!(stap)

    # Update
    sols_x = value.(x)
    sols_y = value.(y)

    iters_x = [sols_x[:]]
    iters_y = [sols_y[:]]
    
    ###
    ### Sequence of penalised sub-problems
    ###
    
    #r_array = [10, 100, 1000, 10000, 100000, 1000000]
    
    for r in r_array
        @objective(stap,
                   Min,
                   p_obj(x, y, r, bounded_indices))
        
        #for multiple ods another loop has to be added
        for k in 1:m
            set_start_value(x[k,1], sols_x[k])
            set_start_value(y[k,1], sols_y[k])
        end           
        
        optimize!(stap)
        
        sols_x = value.(x)
        sols_y = value.(y)

        push!(iters_x, sols_x[:])
        push!(iters_y, sols_y[:])
    end
    
    return iters_x, iters_y
    # Think about returning dual variables
end

"""
Run penalty method for network with restricted flow on some
links for a range of penetration rates
"""
function me_excluded_assignment_penalty_pr(rn, ods, demands, γ_array, bounded_indices, r_array)
    
    xx = []
    yy = []
    agg_a = []
    for γ in γ_array
        x, y = me_excluded_assignment_penalty(rn,
                                              ods,
                                              demands,
                                              γ,
                                              bounded_indices,
                                              r_array)
        push!(xx, x)
        push!(yy, y)
        push!(agg_a, x+y)
    end
    
    return xx, yy, agg_a
end


###
### Def Braess network
###

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

# Optional plotting of networks

#braess net with good path
#ee = [(e.src, e.dst) for e in pe_1]
#save_graph_tikz_edg(rn.g,
#                    [ee],
#                    positions,
#                    "braess_p1.tex",
#                    imp_nodes=[ods[1][1], ods[1][2]],
#                    imp_labels=["O", "D"])
                    
#braess net with bad path
#ee = [(e.src, e.dst) for e in pe_2]
#save_graph_tikz_edg(rn.g,
#                    [ee],
#                    positions,
#                    "braess_p2.tex",
#                    imp_nodes=[ods[1][1], ods[1][2]],
#                    imp_labels=["O", "D"])
                    
#braess net with ugly path
#ee = [(e.src, e.dst) for e in pe_3]
#save_graph_tikz_edg(rn.g,
#                    [ee],
#                    positions,
#                    "braess_p3.tex",
#                    imp_nodes=[ods[1][1], ods[1][2]],
#                    imp_labels=["O", "D"])


# Penalty term array (starting from second iteration)
r_array = [10^i for i in 0:6]

# Actual assignment and optimisation
good_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_1,
                                                 r_array)

bad_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_2,
                                                 r_array)
                                                 
double_results = me_excluded_assignment_penalty_pr(rn,
                                                 ods,
                                                 demands,
                                                 γ_array,
                                                 oi_3,
                                                 r_array)
                                                 
                                                 
#=
The variables good_results and bad_results are arrays of arrays, 
their indices correspond to:
    
    1 - HV flows
    2 - AV flows
    3 - Aggregate flows
    
Each of these has n arrys as elements (one for each value of γ)

Each of these then has one array for each penalty iteration
=#

# Calculating costs

# For γ = 0.5

for γ_ind in [1, 3, 5] 

    # Case 1
    tot_cost = [total_cost(good_results[3][γ_ind][i], a, b) 
                for i in 1:length(r_array)]
                
    hv_cost_per_v = [partial_cost(good_results[1][γ_ind][i],
                                  good_results[3][γ_ind][i],
                                  a,
                                  b) 
                     for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v = [partial_cost(good_results[2][γ_ind][i],
                                  good_results[3][γ_ind][i],
                                  a,
                                  b)
                     for i in 1:length(r_array)] ./ γ_array[γ_ind]*d
                    
    #hv_cost_per_v = [total_cost(good_results[1][γ_ind][i], a, b) 
    #                 for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])
                     
    #av_cost_per_v = [total_cost(good_results[2][γ_ind][i], a, b)
    #                 for i in 1:length(r_array)] ./ γ_array[γ_ind]
    
    # Case 2
    tot_cost2 = [total_cost(bad_results[3][γ_ind][i], a, b) 
                 for i in 1:length(r_array)]
                 
                 
    hv_cost_per_v2 = [partial_cost(bad_results[1][γ_ind][i],
                                  bad_results[3][γ_ind][i],
                                  a,
                                  b) 
                     for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v2 = [partial_cost(bad_results[2][γ_ind][i],
                                  bad_results[3][γ_ind][i],
                                  a,
                                  b)
                     for i in 1:length(r_array)] ./ γ_array[γ_ind]*d
                     
                     #########
                 
    #hv_cost_per_v2 = [total_cost(bad_results[1][γ_ind][i], a, b)
    #                  for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])
                      
    #av_cost_per_v2 = [total_cost(bad_results[2][γ_ind][i], a, b)
    #                  for i in 1:length(r_array)] ./ γ_array[γ_ind]
    
    # Case 3
    tot_cost3 = [total_cost(double_results[3][γ_ind][i], a, b) 
                 for i in 1:length(r_array)]
    
    
    hv_cost_per_v3 = [partial_cost(double_results[1][γ_ind][i],
                                   double_results[3][γ_ind][i],
                                   a,
                                   b) 
                      for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])*d

    av_cost_per_v3 = [partial_cost(double_results[2][γ_ind][i],
                                   double_results[3][γ_ind][i],
                                   a,
                                   b)
                      for i in 1:length(r_array)] ./ γ_array[γ_ind]*d
    
    
    
    #hv_cost_per_v3 = [total_cost(double_results[1][γ_ind][i], a, b)
    #                  for i in 1:length(r_array)] ./ (1 - γ_array[γ_ind])
                      
    #av_cost_per_v3 = [total_cost(double_results[2][γ_ind][i], a, b)
    #                  for i in 1:length(r_array)] ./ γ_array[γ_ind]

    # Plotting

    pen_iters = 0:length(r_array)-1
    
    # Total Costs (all 3 cases)
    
    # Legend labels
    l1 = "Cost (upper protected)"
    l2 = "Cost (interior protected)"
    l3 = "Cost (both outer routes protected)"
    
    plt = plot(pen_iters,
               [tot_cost tot_cost2 tot_cost3],
               label=[l1 l2 l3],
               linewidth = 2,
               markershape=[:utriangle :circle :square],
               markeralpha = 0.6,
               legend=:topleft,
               foreground_color_legend = nothing,
               background_color_legend = nothing,
               color_palette=palette(:tab10))
    
    #ylims!(0.7,0.75)
    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Cost", yguidefontsize=16)
    title!("Total costs (d=$d, γ=$(γ_array[γ_ind]))", titlefontsize=16)
    
    savefig(plt, "total_costs_braess_gamma_d_$(d)_$(γ_array[γ_ind]).pdf")
    
   # Costs per vehicle (cases 2 and 3)
   
   plt2 = plot(pen_iters,
          [hv_cost_per_v2 av_cost_per_v2],
          label=["HVs inner" "AVs inner"],
          linewidth = 2,
          markershape=[:dtriangle :hexagon],
          markeralpha = 0.6,
          legend=:bottomright,
          color_palette=palette(:tab10)[[1,10]])

    #plot!(pen_iters,
    #     [hv_cost_per_v2 av_cost_per_v2],
    #     label=["HVs inner" "AVs inner"],
    #     markershape=[:dtriangle :circle],
    #     legend=:topright,
    #     color_palette=palette(:tab10))
         
    plot!(pen_iters,
         [hv_cost_per_v3 av_cost_per_v3],
         label=["HVs both outer" "AVs both outer"],
         linewidth = 2,
         markershape=[:utriangle :circle],
         markeralpha = 0.6,
         legend=:bottomright,
         foreground_color_legend = nothing,
         background_color_legend = nothing,
         color_palette=palette(:tab10)[[2,4]])
         
    ylabel!("Cost per vehicle", yguidefontsize=16)
    xlabel!("Penalty method iterations", xguidefontsize=16)
    title!("Per-vehicle costs (d=$d, γ=$(γ_array[γ_ind]))", titlefontsize=16)

    savefig(plt2, "perv_costs_braess_gamma_d_$(d)_$(γ_array[γ_ind]).pdf")
end


#γ_ind = 3

#tot_cost = [total_cost(good_results[3][γ_ind][i], a, b) for i in 1:length(r_array)]
#hv_cost_per_v = [total_cost(good_results[1][γ_ind][i], a, b) for i in 1:length(r_array)]
#av_cost_per_v = [total_cost(good_results[2][γ_ind][i], a, b) for i in 1:length(r_array)]

#tot_cost2 = [total_cost(bad_results[3][γ_ind][i], a, b) for i in 1:length(r_array)]
#hv_cost_per_v2 = [total_cost(bad_results[1][γ_ind][i], a, b) for i in 1:length(r_array)]
#av_cost_per_v2 = [total_cost(bad_results[2][γ_ind][i], a, b) for i in 1:length(r_array)]

#tot_cost3 = [total_cost(double_results[3][γ_ind][i], a, b) for i in 1:length(r_array)]
#hv_cost_per_v3 = [total_cost(double_results[1][γ_ind][i], a, b) for i in 1:length(r_array)]
#av_cost_per_v3 = [total_cost(double_results[2][γ_ind][i], a, b) for i in 1:length(r_array)]


# Plot the inner route and both outer routes protected

#pen_iters = 0:length(r_array)-1

#plt2 = plot(pen_iters,
#          [hv_cost_per_v2 av_cost_per_v2],
#          label=["HVs inner" "AVs inner"],
#          markershape=[:dtriangle :circle],
#          legend=:topright,
#          color_palette=palette(:tab10)[[1,10]])

#plot!(pen_iters,
#     [hv_cost_per_v2 av_cost_per_v2],
#     label=["HVs inner" "AVs inner"],
#     markershape=[:dtriangle :circle],
#     legend=:topright,
#     color_palette=palette(:tab10))
     
#plot!(pen_iters,
#     [hv_cost_per_v3 av_cost_per_v3],
#     label=["HVs both outer" "AVs both outer"],
#     markershape=[:dtriangle :circle],
#     legend=:topright,
#     color_palette=palette(:tab10)[[2,4]])
     
#ylabel!("Cost per vehicle")
#xlabel!("Penalty method iterations")
#title!("Per-vehicle costs (d=$d, γ=$(γ_array[γ_ind]))")

#savefig(plt2, "perv_costs_braess_gamma d_$(d)_$(γ_array[γ_ind]).pdf")
