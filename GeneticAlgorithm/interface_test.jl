using NDGeneticAlgorithm,
      TrafficNetworks2,
      SkeletonCities2,
      LightGraphs,
      MetaGraphs,
      GraphIO

###
### Include files for ME-STAP
###

include("../Scratch/tools.jl")
include("../Scratch/penalty.jl")

### Network
g_fname = "test_graph_a_0.5_b_1.5.mg"
rn = load_rn(g_fname)
mg = loadgraph(g_fname, MGFormat())
sorted_edges = collect(edges(rn.g))
num_edges = ne(rn.g)

### Demand
d_range = LinRange(0.0001, 0.03, 10)
γ_range = LinRange(0.0, 1.0, 10)

### OD structure (single OD pair for now)
O = node_closest_to(rn, [0,0]) # lower left-most node
D = node_closest_to(rn, [1,1]) # upper right-most node

ods = [(O,D)]

### MST
mst_inds, mst_edges = find_mst(mg, sorted_edges)
num_mst_edges = length(mst_inds)

### GA params
p_mut = 0.05 # mutation probability

### Genome-edges mapping dictionary
genome_edge_dict = Dict{Int,Int}()
let k = 1
for i in 1:num_edges
    if sorted_edges[i] ∉  mst_edges
        genome_edge_dict[k] =  i
        k += 1
    end
end
end

genome_len = length(genome_edge_dict)


########################################################################
########################################################################


### Construction of initial population
n_ind = 6 # number of individuals in population (even)

pop_ini = Individual[]
for i in 1:n_ind
    push!(pop_ini, Individual(BitArray(rand([0,1], genome_len)), -Inf))
end

ind_1 = pop_ini[1]
fitness!(ind_1, rn, genome_edge_dict, ods, d_range, γ_range)
@show ind_1.fitness

### Calculate fitness and sort
fitness!(pop_ini, rn, genome_edge_dict, ods, d_range, γ_range)
fitness_sort!(pop_ini)

### Select mating pairs (fitness proportional selection)
ma1 = select_mating_array(pop_ini, Int(n_ind/2))
ma2 = select_mating_array(pop_ini, Int(n_ind/2))

# Reproduction
offspring_pop = Individual[]
for k in 1:length(ma1)
    a, b = crossover(ma1[k], ma2[k])
    append!(offspring_pop, [a, b])
end

op1 = Individual[]
for k in 1:length(offspring_pop)
    o_ind = offspring_pop[k]
    new_ind = Individual(copy(o_ind.genome), (o_ind.fitness))
    push!(op1, new_ind)
end

@show op1 .== offspring_pop

# Mutation
for a in offspring_pop
   mutation!(a, p_mut)
end

@show op1 .== offspring_pop
