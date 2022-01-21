### Type for Individuals in the population for GA


# Individuals (maybe individuals should be immutable...)

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

"""
Generate a random population of `m`  individuals with 
length of genome: `genome_len`.
"""
function generate_population(genome_len, m)
    pop = Individual[]
    for i in 1:m
        push!(pop, Individual(BitArray(rand([0,1], genome_len)), -Inf))
    end
    pop
end


### Function for mapping genome to non-mst edges

function genome_edge_mapping(sorted_edges, mst_edges)
    num_edges = length(sorted_edges)
    genome_edge_dict = Dict{Int,Int}()

    k = 1
    for i in 1:num_edges
        if sorted_edges[i] ∉  mst_edges
            genome_edge_dict[k] =  i
            k += 1
        end
    end
    genome_edge_dict
end
