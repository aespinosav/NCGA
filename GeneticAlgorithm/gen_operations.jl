#=
using Distributions
=#


### Naive genetic operations

"""
Carries out single-point crossover.
If crossover point is provided, do deterministically.
"""
function crossover(a::T, b::T; x_point=nothing) where {T<: AbstractIndividual}

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

    return T(ng1), T(ng2)
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
function mutate!(a::T, mutation_vect::BitArray) where {T<:AbstractIndividual}
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
function mutation!(a::T, p) where {T<:AbstractIndividual}
    N = length(a.genome)
    mv = draw_bernoulli_vec(N, p)
    mutate!(a, mv)
end
