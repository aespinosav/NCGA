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

"""
Generate individual with genome length `gen_len`
and average density of `1`s in genome of `p`
"""
function Individual(gen_len::Int, p::Float64)
    gen_dist = Binomial(1, p)
    genome = BitArray(rand(gen_dist, gen_len))
    Individual(genome)
end


function ==(a::Individual, b::Individual)
    a.genome == b.genome &&
    a.fitness == b.fitness
end

length(a::Individual) = length(a.genome)

# Population: for now to be treated as an array/set of individuals

"""
Generate a random population of `m`  individuals with 
length of genome: `genome_len`.

ρ: Average density of `1` in genome
"""
function generate_population(genome_len, m; ρ=0.5)
    pop = Individual[]
    for i in 1:m
        a = Individual(genome_len, ρ)
        push!(pop, a)
    end
    pop
end

"""
Generates population with "genetic" diversity
by drawing genomes with different densities of `1`s
scanning the densiity range `dens_range`
"""
function generate_diverse_pop(genome_len, pop_n; dens_min=0.0, dens_max=1.0)
    densities = LinRange(dens_min, dens_max, pop_n)

    pop = Individual[]
    for d in densities
        a = Individual(genome_len, d)
        push!(pop, a)
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
