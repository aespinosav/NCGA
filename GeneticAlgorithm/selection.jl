# Functions for selecting mating pairs

#=
using StatsBase
=#

### Fitness proportionate selection
"""
Returns array with probability of selection (proportional to fitness)
in the order individuals are in the population array
"""
function fitprop_selection_prob(pop::Array{Individual,1})
    fitnesses = [a.fitness for a in pop]
    p_selection = fitnesses / sum(fitnesses)
    Weights(p_selection)
end

"""
Select single individual randomly proportional to fitness.

This is inefficient since it calculates probabilities of 
all individiuals in pop to select only 1...
"""
function select_individual(pop::Array{Individual,1})
    weights = fitprop_selection_prob(pop)
    sample(pop, weights)
end

function select_individual(pop::Array{Individual,1}, w::AbstractWeights)
    sample(pop, w)
end


"""
Select (returns array) `k` individuals, selected randomly
proportional to fitness
"""
function select_mating_array(pop::Array{Individual,1}, k::Int)
    weights = fitprop_selection_prob(pop)
    sample(pop, weights, k)
end

###
### TODO: alternative selection methods
###

### Tournament (ordinal selection)

"""
Select mating individual from a `q`-tournament
"""
function tournament_selection(pop::Array{Individual,1}, q)
    contestants = rand(pop, q)
    sort!(contestants, by=x->x.fitness, rev=true)
    contestants[1]
end

"""
Select (returns array) `k` individuals, selected via tournamets
of `q` participants.
"""
function select_mating_array(pop::Array{Individual,1}, k::Int, q::Int)
    mating_array = Individual[]
    for i in 1:k
        victor = tournament_selection(pop, q)
        push!(mating_array, victor)
    end
    mating_array
end
