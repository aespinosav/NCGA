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
       Base.length,
       TrafficNetworks2.multi_pair_stap_nc

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
       dispersion_fitness!,
       fitness_sort!,
       stap_wrapper_fit,
       # From selection.jl
       select_mating_array,
       # From evolve.jl
       evolve,
       elite_conservation,
       introduce_diverse_individuals!,
       # From plotting.jl
       plot_fit_dens,
       plot_evolution,
       term_link_plot,
       # From stap_tools.jl
       multi_pair_stap_nc,
       node_closest_to,
       find_mst,
       av_links,
       travel_times,
       marginal_travel_times,
       total_cost,
       partial_cost

# Module files

include("chromo.jl")
include("gen_operations.jl")
include("fitness.jl")
include("selection.jl")
include("evolve.jl")
include("plotting.jl")
include("stap_tools.jl")
end
