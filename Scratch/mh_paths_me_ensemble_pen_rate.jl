#=
Script for plotting results for ensembles as functions
of penetration rate and demand
=#

using Plots,
      FileIO,
      DelimitedFiles

###
### Main
###

function main()

    ###
    ### Read parameter file
    ###

    params = load("batch_run_params.jld2")
    # manually for now (should save them in param file!)
    α_hat = 0.75
    β = 1.5
    # add mean number of restricted links?

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

    # Mean UE cost
    mean_ue_cost = params["ue_cost"]


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

    # Normalised costs
    norm_mean_tot_costs = readdlm("normed_mean_total_costs.dat")
    norm_std_tot_costs = readdlm("std_normed_total_costs.dat")

    # Per-vehicle cost differences
    pv_cost_diff = readdlm("perv_mean_costs_diff.dat")
    std_pv_cost_diff =  readdlm("perv_std_cost_diff.dat")

    ### Select results after penalty convergence

    tc =  mean_tot_costs[end,:]
    s_tc = std_tot_costs[end,:]
    ntc = norm_mean_tot_costs[end,:]
    s_ntc = norm_std_tot_costs[end,:]
    mpc_hv = mean_perv_costs_hv[end,:]
    s_mpc_hv = std_perv_costs_hv[end,:]
    mpc_av = mean_perv_costs_av[end,:]
    s_mpc_av = std_perv_costs_av[end,:]

    #pvc_diff = mpc_hv - mpc_av
    #pvc_diff ./= mpc_hv # how many times the hv cost...
    pvc_diff = mpc_hv ./ mpc_av

    # Plot tot cost
    plt1 = plot(γ_array,
                tc,
                ribbon=s_tc,
                label="System cost",
                markershape=:circle)
    xlims!(0,1)
    ylabel!("Mean cost", yguidefontsize=16)
    xlabel!("Penetration rate, γ", xguidefontsize=16)
    title!("$(nn) Networks, d=$d, α̂=$(α_hat), β=$(β)", titlefontsize=14)
    savefig(plt1, "total_cost_with_gamma.pdf")


    # Plot normed cost
    plt2 = plot(γ_array,
                ntc,
                ribbon=s_ntc,
                label="Normalised system cost",
                markershape=:circle,
                legend=false)
    hline!([1], lc=:grey, linestyle=:dash)
    xlims!(0,1)
    lens!([0.5, 1], [0, 1], inset = (1, bbox(0.5, 0.0, 0.4, 0.4)), subplot=2)
    ylabel!("Mean normalised cost", yguidefontsize=16, subplot=1)
    xlabel!("Penetration rate, γ", xguidefontsize=16, subplot=1)
    title!("$(nn) Networks, d=$d, α̂=$(α_hat), β=$(β)", titlefontsize=14, subplot=1)
    savefig(plt2, "norm_cost_with_gamma.pdf")


    # Plot per-veh cost
    plt3 = plot(γ_array,
               [mpc_hv mpc_av],
               ribbon=[s_mpc_hv s_mpc_av],
               label=["HV pv cost" "AV pv cost"],
               markershape=[:dtriangle :circle],
               markeralpha = 0.6,
               color_palette=palette(:tab10)[[1,10]])
    hline!([mean_ue_cost/d], lc=:red, linestyle=:dash, label="UE pv cost")
    xlims!(0,1)
    ylabel!("Per-vehicle cost", yguidefontsize=16)
    xlabel!("Penetration rate, γ", xguidefontsize=16)
    title!("$(nn) Networks, d=$d, α̂=$(α_hat), β=$(β)", titlefontsize=14)
    savefig(plt3, "pv_cost_with_gamma.pdf")


    # Per-veh cost difference
    plt4 = plot(γ_array,
                pvc_diff,
                label="Per-veh cost difference (HV - AV)",
                markershape=[:circle],
                legend=:bottomright)
    hline!([0], linecolor=:grey, linestyle=:dash)
    xlims!(0,1)
    ylabel!("Per-vehicle cost difference", yguidefontsize=16)
    xlabel!("Penetration rate, γ", xguidefontsize=16)
    title!("$(nn) Networks, d=$d, α̂=$(α_hat), β=$(β)", titlefontsize=14)
    savefig(plt4, "pv_cost_diff_with_gamma.pdf")

end

main()
