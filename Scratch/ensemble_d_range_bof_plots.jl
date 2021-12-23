#=
This script plots the ensemble averaged costs of the
network control scheme calculated in `me_ensemble_d_range_bof.jl`
and saved in the corresponding data files
=#

using Plots,
      FileIO

function main()

    ###
    ### Load data and range defs
    ###

    # Will not define vars explicitly (use the dict!)
    data = load("d_range_bof_data.jld2")

    γ_range = data["gamma_init"]:data["gamma_step"]:data["gamma_stop"]
    d_range = data["d_start"]:data["d_step"]:data["d_stop"]

    ###
    ### Plots
    ###

    ### Total cost plot
    labs = reshape(["γ = $(γ)" for γ in γ_range],
                    1,
                    length(γ_range))

    plt1 = plot(d_range,
                data["mean_tot_costs"],
                #ribbon = data["std_tot_costs"],
                label = labs,
                markershape = [:circle :dtriangle :square],
                legend=:topleft)

    plot!(d_range,
          [data["mean_ue"] data["mean_so"]],
          #ribbon=[data["std_ue"] data["std_so"]],
          label = ["UE" "SO"],
          linestyle = :dash,
          linecolor = [:grey :black])

    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Mean ensemble cost", yguidefontsize=16)

    savefig(plt1, "tot_ens_costs_d_range_bof.pdf")

    ### Normalised cost plot

    max_dat = maximum(data["mean_norm_costs"])

    plt2 = plot(d_range,
                data["mean_norm_costs"],
                ribbon = data["std_norm_costs"],
                label = labs,
                markershape = [:circle :dtriangle :square])

    hline!([1],
           label="UE",
           lc=:grey,
           linestyle=:dash)

    hline!([0],
            label="SO",
            linecolor=:black,
            linestyle=:dash)

    ylims!(-0.01 ,maximum([1.2, 1.1*max_dat]))
    xlims!(d_range[2], 0.02) # a mano el top lim...

    xlabel!("Demand", xguidefontsize=16)
    ylabel!("Mean normalised ensemble cost", yguidefontsize=16)

    savefig(plt2, "norm_ens_costs_d_range_bof.pdf")

    ### Per-vehicle cost plot
end

main()
