#=
This script plots the ensemble averaged costs of the
network control scheme calculated in `me_ensemble_d_range_bof.jl`
and saved in the corresponding data files
=#

using Plots,
      FileIO,
      DelimitedFiles


function main()

    # Will not define vars explicitly (use the dict!)
    data = load("d_range_bof_data.jld2")

    γ_range = data["γ_init"]: data["γ_step"]: data["γ_stop"]
    d_range = data["d_start"]:data["d_step"]:data["d_stop"]

    


end
