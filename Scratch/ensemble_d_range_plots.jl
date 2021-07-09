#=
This script plots the ensemble averaged costs of the
network control scheme calculated in `me_ensemble_d_range.jl`
and saved in the corresponding data files
=#

using Plots,
      FileIO,
      DelimitedFiles

function main()

    ###
    ### Read parameter file
    ###

    params = load("batch_run_params_d_range.jld2")

    # manually for now (should save them in param file!)
    α_hat = 0.75
    β = 1.5

    # Reconstruct paramters
    optimiser = params["optim"]

    # Penetration rate range
    γ_init = params["gamma_init"]
    γ_step = params["gamma_step"]
    γ_stop = params["gamma_stop"]
    γ_array = γ_init:γ_step:γ_stop

    # Demand range
    d_start = params["d_start"]
    d_step = params["d_step"]
    d_stop = params["d_stop"]
    demand_range = d_start:d_step:d_stop

    # MH parameters
    μ = params["μ"]
    p_splice = params["p_splice"]

    # Dirs and net files
    ens_dir = params["dir"]
    nn = params["num_nets"]

    ###
    ### Read data files
    ###

    ue_costs = readdlm("d_range_mean_ue_costs.dat")
    so_costs = readdlm("d_range_mean_so_costs.dat")

    tot_costs_mean = readdlm("d_range_mean_tot_costs.dat")
    tot_costs_std = readdlm("d_range_std_tot_costs.dat")

    normed_costs_mean = readdlm("d_range_mean_normed_tot_costs.dat")
    normed_costs_std = readdlm("d_range_std_normed_tot_costs.dat")

    per_veh_hv_mean = readdlm("d_range_mean_perv_costs_hv.dat")
    per_veh_hv_std = readdlm("d_range_std_perv_costs_hv.dat")

    per_veh_av_mean = readdlm("d_range_mean_perv_costs_av.dat")
    per_veh_av_std = readdlm("d_range_std_perv_costs_av.dat")


    ###
    ### Plot
    ###

    ### Total cost plot

    labs = reshape(["γ = $(γ)" for γ in γ_array],
                    1,
                    length(γ_array))

    plt1 = plot(demand_range,
                tot_costs_mean,
                ribbon = tot_costs_std,
                label = labs,
                markershape = [:circle :dtriangle :square],
                legend=:topleft)

    plot!(demand_range,
          [ue_costs[:,1] so_costs[:,1]],
          #ribbon=[ue_costs[:,2] so_costs[:,2]],
          label = ["UE costs" "SO costs"],
          linestyle = :dash,
          linecolor = [:grey :black])

    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Mean ensemble cost", yguidefontsize=16)

    savefig(plt1, "tot_ens_costs_d_range.pdf")


    ### Normed cost plot

    plt2 = plot(demand_range,
                normed_costs_mean,
                #ribbon = normed_costs_std,
                label = labs,
                markershape = [:circle :dtriangle :square])

    hline!([1, 0],
           label=["UE cost" "SO cost"],
           lc=[:grey :black],
           linestyle=:dash)

    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Mean normalised ensemble cost", yguidefontsize=16)

    savefig(plt2, "norm_ens_costs_d_range.pdf")

    # Per-veh cost plot
end

main()
