# Initial testing for Genetic Algorithm for network design

module NDGeneticAlgorithm

using TrafficNetworks2,
      SkeletonCities2,
      Graphs,
      GraphIO,
      MetaGraphs,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Ipopt,
      StatsBase,
      Distributions,
      UnicodePlots

import Base.==,
       Base.length

export # From chromo.jl
       Individual,
       length,
       genome_edge_mapping,
       generate_population,
       generate_diverse_pop,
       # From gen_operations.jl
       crossover,
       mutation!,
       # From fitness.jl
       fitness!,
       fitness_sort!,
       # From selection.jl
       select_mating_array,
       # From evolve.jl
       evolve,
       elite_conservation,
       # From plotting.jl
       plot_fit_dens,
       plot_evolution

# Module files

include("chromo.jl")
include("gen_operations.jl")
include("fitness.jl")
include("selection.jl")
include("evolve.jl")

end
