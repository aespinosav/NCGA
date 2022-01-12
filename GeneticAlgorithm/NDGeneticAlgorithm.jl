# Initial testing for Genetic Algorithm for network design


module NDGeneticAlgorithm

using Distributions

import Base.==,
       Base.length

export Individual,
       length,
       fitness,
       crossover,
       mutation!

###
### Types and associated functions
###

# Individuals

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
Individual(G::BitArray) = Individual(G, -Inf)

function ==(a::Individual, b::Individual)
    a.genome == b.genome &&
    a.fitness == b.fitness
end

length(a::Individual) = length(a.genome)


# Population: for now to be treated as an array/set of individuals


###
### Functions
###

"""
Calculates fitness of individual
"""
function fitness(a::Individual, params)

end


### Naive genetic operations

"""
Carries out single-point crossover.
If crossover point is provided, do deterministically.
"""
function crossover(a::Individual, b::Individual; x_point=nothing)

    A = a.genome
    B = b.genome

    len = length(A)

    ng1 = BitArray(undef, len)
    ng2 = BitArray(undef, len)

    if x_point == nothing
        x_point = rand(2:len-1)
    end

    ng1[1:x_point] = A[1:x_point]
    ng1[x_point:end] = B[x_point:end]

    ng2[1:x_point] = B[1:x_point]
    ng2[x_point:end] = A[x_point:end]

    return Individual(ng1), Individual(ng2)
end

"""
Draws random bernoulli variable vector
for mutate! function
"""
function draw_bernoulli_vec(N::Int, p::Float64)
    dist = Binomial(1, p)
    vec = BitArray(rand(dist, N))
end

"""
Mutate genome according to `mutation vector`
"""
function mutate!(a::Individual, mutation_vect::BitArray)
    N = length(a.genome)
    for i in 1:N
        if mutation_vect[i] == 1
            # Filp bit in genome
            a.genome[i] = !(a.genome[i])
        end
    end
end

"""
Function that mutates an individual's (`a`) genome
according to probability `p`.

Calls mutate! after generating the mutation vector
and passing it as an argument.
"""
function mutation!(a::Individual, p)
    N = length(a.genome)
    mv = draw_bernoulli_vec(N, p)
    mutate!(a, mv)
end


# Module end
end
