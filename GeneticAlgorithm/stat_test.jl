#=
using StatsBase,
      Distributions,
      HypothesisTests
=#

"""
Performs a two sample equal variance t-test on two samples
(of size `n0`) from an array of 'random variables' `var_array`
with confidence level `alpha`.

Usage
------

    test_same_means(var_array, n0, alpha)

Returns
-------
    + `true`: if sample means are "the same" (accept H₀: X̄₁=X̄₂)
    + `false`: false if they are not (reject H₀: X̄₁≠X̄₂)
"""
function test_same_means(var_array, n0, alpha)

    # Index range
    r = 1:length(var_array)

    #Take samples
    sample1 = sample(r, n0, replace=false)
    rem = findall(x->!(x in sample1), r)
    sample2 = sample(rem, n0, replace=false)

    X₁ = var_array[sample1]
    X₂ = var_array[sample2]

    # Test
    t_test = EqualVarianceTTest(X₁, X₂)
    p_value = pvalue(t_test, tail=:both)

    p_value > alpha
end
