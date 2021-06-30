#=
    Script for plotting output of `mh_paths_me_ensemble.jl`
=#

using Plots,
      FileIO,
      DelimitedFiles


###
### Read parameter file
###

params = load("batch_run_params.jld2")

# Reconstruct paramters
optimiser = params["optim"]

# Penalty parameters
penalty_steps = params["penalty_steps"]
r_array = [10^i for i in 1:penalty_steps]

γ_init = params["gamma_init"]
γ_step = params["gamma_step"]
γ_stop = params["gamma_stop"]
# Penetration rate 
γ_array = γ_init:γ_step:γ_stop

# MH parameters
μ = params["μ"]
d = params["d"]
p_splice = params["p_splice"]

# Dirs and net files
ens_dir = params["dir"]
nn = params["num_nets"]


###
### Read data files
###

mean_tot_costs = readdlm("mean_total_costs.dat")
std_tot_costs = readdlm("std_total_costs.dat")

# HVs
mean_perv_costs_hv = readdlm("mean_perv_costs_hv.dat")
std_perv_costs_hv = readdlm("std_perv_costs_hv.dat")
# AVs
mean_perv_costs_av = readdlm("mean_perv_costs_av.dat")
std_perv_costs_av = readdlm("std_perv_costs_av.dat")

norm_mean_tot_costs = readdlm("normed_mean_total_costs.dat")
norm_std_tot_costs = readdlm("std_normed_total_costs.dat")

###
### Plotting
###

pen_iters = 0:length(r_array)

y_min, y_max = minimum(mean_tot_costs), maximum(mean_tot_costs)
std_margin = maximum(std_tot_costs)
data_margin = 0.15*(y_max - y_min)
margin = maximum([std_margin, data_margin])

ny_min, ny_max = 0, maximum(norm_mean_tot_costs)
margin2 = 0.1*ny_max             

# Total cost plot

for (i,γ) in enumerate(γ_array)

    # Total cost plot
    plt = plot(pen_iters,
               mean_tot_costs[:,i],
               ribbon=std_tot_costs[:,i],
               label="Total cost",
               markershape=[:circle],
               markeralpha = 0.6,
               legend=:bottomright,
               color_palette=palette(:tab10))

    ylims!(y_min - margin, y_max + margin)

    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Cost", yguidefontsize=16)
    title!("Total costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt, "ensemble_mean_total_costs_gamma_$(γ).pdf")

    ##########################################################

    # Per-vehicle cost plot
    plt2 = plot(pen_iters,
               [mean_perv_costs_hv[:,i] mean_perv_costs_av[:,i]],
               ribbon=[std_perv_costs_hv std_perv_costs_av],
               label=["HV pv cost" "AV pv cost"],
               markershape=[:dtriangle :circle],
               markeralpha = 0.6,
               color_palette=palette(:tab10)[[1,10]])

    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Class cost (per vehicle)", yguidefontsize=16)
    title!("Per-vehicle costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt2, "ensemble_mean_pv_costs_gamma_$(γ).pdf")

    ###########################################################

    # Normalised costs
    plt3 = plot(pen_iters,
                norm_mean_tot_costs[:,i],
                ribbon=norm_std_tot_costs,
                label="Normalised tot. cost",
                markershape=[:circle],
                markeralpha = 0.6,
                legend=:bottomright,
                color_palette=palette(:tab10))

    ylims!(ny_min, ny_max + margin2)

    xlabel!("Penalty method iterations", xguidefontsize=16)
    ylabel!("Normalised cost", yguidefontsize=16)
    title!("Normalised total costs (γ=$(γ), d=$(d))", titlefontsize=16)
    savefig(plt3, "ensemble_normed_mean_total_costs_gamma_$(γ).pdf")

    # Cost difference

end
