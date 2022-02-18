# Functions for evolving the population from generation to generation

"""
Select mating pairs from `pop` and produce generation of offspring.
Includes mutation step with each allele having probability
`p_mut` of flipping.

The `select_mating_array` method depends on the arguments
used to call it.

If `q`== `nothing` then fitness proportional selection is
used, otherwise an `Int` should be passed that is the number
of individuals per tournament for `tournament_selection`
"""
function evolve(pop::Array{Individual,1}, p_mut; q=nothing)

    n_ind = length(pop)
    n_mat_pairs = round(Int, n_ind/2)

    if q==nothing
        ### Select mating pairs (fitness proportional selection)
        ma1 = select_mating_array(pop, n_mat_pairs)
        ma2 = select_mating_array(pop, n_mat_pairs)
    elseif typeof(q) <: Int
        ma1 = select_mating_array(pop, n_mat_pairs, q)
        ma2 = select_mating_array(pop, n_mat_pairs, q)
    end

    # Reproduction
    offspring_pop = Individual[]
    for k in 1:length(ma1)
        a, b = crossover(ma1[k], ma2[k])
        append!(offspring_pop, [a, b])
    end

    # Mutation
    for a in offspring_pop
       mutation!(a, p_mut)
    end

    offspring_pop
end

"""
Keep top `frac` prcent of parent population in the next generation
as long as they are fitter than the fittest in offspring generation

Assumes that pop is sorted by fitness.
"""
function elite_conservation(pop, off_pop, frac)

    n_pop = length(pop)
    ind_frac = round(Int, frac*n_pop)

    # Possible elites
    elite_pop = pop[1:ind_frac]

    top_off = off_pop[1].fitness

    cutoff = 0
    for (k, ind) in enumerate(elite_pop)
        if ind.fitness > top_off
            cutoff = k
        end
    end

    if cutoff > 0
        new_pop = [elite_pop[1:cutoff]..., off_pop[1:end-cutoff]...]
    else
        new_pop = off_pop
    end

    new_pop
end


"""
Introduce diverse individuals instead of the 
lowest fitness current individual in the population.

Assume population is already sorted by fitness.
"""
function introduce_diverse_individuals!(pop::Array{Individual,1},
                                        replacement_ratio,
                                        rn,
                                        genome_link_dict,
                                        ods,
                                        d_range_or_float,
                                        γ_range_or_float)

    pop_size = length(pop)
    genome_length = length(pop[1].genome)
    num_new_inds = round(Int, pop_size*replacement_ratio)

    densities = [sum(a.genome)/genome_length for a in pop]
    maximum_density = maximum(densities)

    new_inds = generate_diverse_pop(genome_length,
                                    num_new_inds,
                                    dens_min = maximum_density,
                                    dens_max = 1.0)

    fitness!(new_inds,
             rn,
             genome_link_dict,
             ods,
             d_range_or_float,
             γ_range_or_float)

    pop[end-(num_new_inds-1):end] = new_inds
end
