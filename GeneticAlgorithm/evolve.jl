# Functions for evolving the population from generation to generation

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

