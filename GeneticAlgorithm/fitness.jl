# Functions for calculating fitness and interfacing with STAP code

#=
using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      GraphIO,
      MetaGraphs,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Ipopt
=#

"""
Calculates fitness of individual.

This fitness function does a very coarse integration of the total costs
in (d, γ)-space for the given demand and penetration rate ranges.

Very naive method, but we can see if it works.
"""
function fitness!(a::Individual,
                 rn::RoadNetwork,
                 genome_link_dict,
                 ods,
                 d_range,
                 γ_range)

    n_ods = length(ods)

    ld = length(d_range)
    lg = length(γ_range)
    num_gridpoints = ld * lg

    cost_array = zeros(ld, lg)

    for (i,q) in enumerate(d_range)
        for (j,γ) in enumerate(γ_range)

            demands = (q/n_ods) * ones(n_ods)

            true_alleles = findall(x->x==1, a.genome)
            bound_links = [genome_link_dict[k] for k in true_alleles]

            x, y = stap_wrapper_fit(rn,
                                    ods,
                                    demands,
                                    γ,
                                    bound_links)

            agg = x + y

            tc = total_cost(agg, rn.edge_params[:a], rn.edge_params[:b])
            tc_av = tc / q

            cost_array[i,j] = tc_av
        end
    end
    fit = 1 / (sum(cost_array) / num_gridpoints)
    a.fitness = fit
end

function fitness!(pop::Array{Individual,1},
                  rn::RoadNetwork,
                  genome_link_dict,
                  ods,
                  d_range,
                  γ_range)

    for ind in pop
        fitness!(ind, rn, genome_link_dict, ods, d_range, γ_range)
    end
end


"""
Sorts a population based on fitness ranking.
Descending order
"""
function fitness_sort!(pop::Array{Individual,1})
    sort!(pop, by= p -> p.fitness, rev=true)
end


"""
ods = [(O₁, D₁), ..., (Oₙ, Dₙ)]
demands: vector of demands for each OD pair
"""
function stap_wrapper_fit(rn, ods, demands, γ, bounded_indices; solver=:ipopt)

    # Structure
    n = nv(rn.g)
    m = ne(rn.g)
    A = incidence_matrix(rn.g)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])

    n_ods = length(ods)

    # Vectors for conservation constraints
    d_vects_hv = SparseVector{Float64,Int64}[]
    d_vects_av = SparseVector{Float64,Int64}[]
    for i in 1:length(ods)
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]

        push!(d_vects_hv, (1-γ).*d_vec)
        push!(d_vects_av, γ.*d_vec)
    end

    # Objective (Mixed multi-class Beckmann UE/SO)
    me_obj_func(x, y) = sum( B[i,i]*(x[i]+y[i])^2 +
                             a[i]*x[i] +
                             0.5*a[i]*y[i] for i in 1:length(x))

    # Choose solver
    if solver == :gurobi
        stap = Model(() -> Gurobi.Optimizer(env))
    elseif solver == :ipopt
        stap = Model(Ipopt.Optimizer)
    end

    # Supress output from solver
    set_silent(stap)

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)

    # OD specific link flows AV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    # Conservation constraints
    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])

    #Bound constraints (HV flow set to 0 on some links)
    for i in bounded_indices
        for k in 1:n_ods
            @constraint(stap, x[i,k]==0)
        end
    end


    # Unpenalised objective
    @objective(stap, Min, me_obj_func(x, y))

    # Solve
    optimize!(stap)
    return value.(x)[:], value.(y)[:]
end
