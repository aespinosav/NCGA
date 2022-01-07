# Initial testing for Genetic Algorithm for network design


module NDGeneticAlgorithm

export Individual, fitness, crossover, mutate


###
### Definitions
###

"""
Individual in the genome population.

`genome`: Array of binary variables (AV exclusive links)
`fitness`: Total system cost (upper level objective)
"""
mutable struct Individual
    genome::BitArray
    fitness::Float64
end
Individual(N::Int) = Individual(BitArray(undef, N), -Inf)


"""
Calculates fitness of individual
"""
function fitness(indi::Individual, params)

end

###
### Naive genetic operations
###

"""
Carries out crossover
"""
function crossover()

end

"""
Mutate offspring genome
"""
function mutate()

end

end
