# Functions for evolving the population from generation to generation

"""
Select mating pairs from `pop` and produce generation of offspring.
Includes mutation step with each allele having probability
`p_mut` of flipping.
"""
function evolve(pop::Array{Individual,1}, p_mut)

    n_ind = length(pop)
    n_mat_pairs = round(Int, n_ind/2)

    ### Select mating pairs (fitness proportional selection)
    ma1 = select_mating_array(pop, n_mat_pairs)
    ma2 = select_mating_array(pop, n_mat_pairs)

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
