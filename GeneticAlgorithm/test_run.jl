using TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO,
      NDGeneticAlgorithm

### GA params

p_mut = 0.05 # mutation probability
n_pop = 50 # Individuals in population (at any generation)
n_gen = 40 # Number of generations


#######################################################
######################################################


###
### Include files for ME-STAP
###

include("../Scratch/tools.jl")
include("../Scratch/penalty.jl")


###
### Network
###

g_fname = "test_graph_a_0.5_b_1.5.mg"
rn = load_rn(g_fname)
mg = loadgraph(g_fname, MGFormat())
sorted_edges = collect(edges(rn.g))
num_edges = ne(rn.g)

### OD structure (single OD pair for now)

O = node_closest_to(rn, [0,0]) # lower left-most node
D = node_closest_to(rn, [1,1]) # upper right-most node

ods = [(O,D)]

### MST

mst_inds, mst_edges = find_mst(mg, sorted_edges)
num_mst_edges = length(mst_inds)

### Demand

d_range = LinRange(0.0001, 0.03, 10)
γ_range = LinRange(0.0, 1.0, 10)


###
### GA setup
###

genome_edge_dict = genome_edge_mapping(sorted_edges, mst_edges)
genome_len = length(genome_edge_dict)

### Starting population

pop = generate_population(genome_len, n_pop)

# Calculate fitness and sort
fitness!(pop, rn, genome_edge_dict, ods, d_range, γ_range)
fitness_sort!(pop)

###
### Evolve population
###

gens = [pop]
for i in 1:n_gen

    # Offspring
    off_pop = evolve(gens[i], p_mut)

    fitness!(off_pop, rn, genome_edge_dict, ods, d_range, γ_range)
    fitness_sort!(off_pop)

    # Elitism (keep top 10% of individuals from previous gen if fitter)
    elite_pop = pop[1:round(Int, 0.1*n_pop)]

    cutoff = 0
    for (k, ind) in enumerate(elite_pop)
        if ind.fitness > off_pop[1].fitness
            cutoff = k
        end
    end

    if cutoff > 0
        new_pop = [elite_pop[1:cutoff]..., off_pop[cutoff:end]...]
    else
        new_pop = off_pop
    end

    # New generation
    push!(gens, new_pop)
end
