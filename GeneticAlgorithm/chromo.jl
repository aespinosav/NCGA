### Type for Individuals in the population for GA


"""
Abstract type for individuals of paopulations used
for the genetic algorithm
"""
abstract type AbstractIndividual end

### Concrete individual types

"""
Individual in the genome population.

`genome`: Array of binary variables (AV exclusive links)
`fitness`: Total system cost (upper level objective)
"""
mutable struct Individual <: AbstractIndividual
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

# Specialised Individual for KSP (has fields for more details on fitness)

"""
A more specialized class of Individual that can store details
about the paths used to caluclate the fitness via the
`dispersion_fitness_ksp"` function.
"""
mutable struct Individual_KSP <: AbstractIndividual
    genome::BitArray
    fitness::Float64
    # For use storing relevant things during fitness calcs
    mean_path_cost::Float64
    min_path_cost::Float64
    max_path_cost::Float64
    num_k_paths::Int
end
Individual_KSP(N::Int) = Individual_KSP(BitArray(undef, N), -Inf, -Inf, -Inf, -Inf, -1)
Individual_KSP(G::BitArray) = Individual_KSP(G, -Inf, -Inf, -Inf, -Inf, -1)

function Individual_KSP(gen_len::Int, p::Float64)
    gen_dist = Binomial(1, p)
    genome = BitArray(rand(gen_dist, gen_len))
    Individual_KSP(genome)
end


###
### Methods (interface) for AbstractIndividual
###

# Requires import from Base in module file

function ==(a::T, b::T) where {T<:AbstractIndividual}
    a.genome == b.genome &&
    a.fitness == b.fitness
end

length(a::T) where {T<:AbstractIndividual} = length(a.genome)


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


###
### Other functions
###

### Function for mapping genome to non-mst edges

"""
Returns a dictionary that maps edges in the genome of individuals
to edges on the network that can be selected as AV exclusive netwroks
(i.e. edges not in the mst)

The dictionary uses the indices of the edges, but takes
in arrays of Edge objects

Usage
----
    genome_edge_mapping(sorted_edges, mst_edges)
"""
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


"""
Generates population with "genetic" diversity
by drawing genomes with different densities of `1`s
scanning the densiity range `dens_range`

I have comandeered it to generate a population of
Individual_KSP

The function should probably take a type as argument
and return population of that type...
"""
function generate_diverse_pop(genome_len, pop_n; dens_min=0.0, dens_max=1.0)
    densities = LinRange(dens_min, dens_max, pop_n)

    pop = Individual_KSP[]
    for d in densities
        a = eltype(pop)(genome_len, d)
        push!(pop, a)
    end
    pop
end

#=
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
=#



