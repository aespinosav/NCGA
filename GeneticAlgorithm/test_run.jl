using TrafficNetworks2,
      SkeletonCities2,
      Graphs,
      MetaGraphs,
      GraphIO,
      Serialization,
      Dates,
      NDGeneticAlgorithm

### GA params

p_mut = 0.05 # mutation probability
n_pop = 10 # Individuals in population (at any generation)
n_gen = 5 # Number of generations

elite_prop = 0.15

### Saving

time_of_run = Dates.format(now(), dateformat"yyyy_mm_dd_HH_MM")

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

# Genome-Net mapping
genome_edge_dict = genome_edge_mapping(sorted_edges, mst_edges)
genome_len = length(genome_edge_dict)

# Starting population
#pop = generate_population(genome_len, n_pop)

pop = generate_diverse_pop(genome_len, n_pop)

# Seeding (this is done automatically by using
# generate_diverse_pop function)

#seed_ind_1 = Individual(BitArray(ones(genome_len)), -Inf)
#seed_ind_2 = Individual(BitArray(zeros(genome_len)), -Inf)
#pop[1] = seed_ind_1
#pop[2] = seed_ind_2

# Calculate fitness and sort
fitness!(pop, rn, genome_edge_dict, ods, d_range, γ_range)
fitness_sort!(pop)


###
### Evolve population
###

gens = [pop]
for i in 1:n_gen

    # Offspring
    #off_pop = evolve(gens[i], p_mut)
    off_pop = evolve(gens[i], p_mut, q=2) # tournament selection

    # Elitism
    new_pop = elite_conservation(pop, off_pop, elite_prop)

    # Fitness calculation
    fitness!(off_pop, rn, genome_edge_dict, ods, d_range, γ_range)
    fitness_sort!(off_pop)

    # New generation
    push!(gens, new_pop)
end

### Save population
serialize("pop_savefile_$(time_of_run).jlsl", gens)
# Use deserialize to open
