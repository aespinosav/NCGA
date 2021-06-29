using SkeletonCities2, LightGraphs, MetaGraphs

n = 8
α_hat = 0.8
β = 1.5

α_crit = 1/(n+1)
α = α_hat*α_crit

net = αβ_network_meta(n, α, β)

savegraph("g_test_medium.mg", net)
