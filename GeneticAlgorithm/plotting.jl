#=
using UnicodePlots
=#

function plot_fit_dens(pop)
    g_len = length(pop[1].genome)
    fitss = [a.fitness for a in pop]
    dens = [sum(a.genome)/g_len for a in pop]

    xl = (0,1)
    yl = (0, 1.15*maximum(fitss))

    scatterplot(dens, fitss, xlim=xl, ylim=yl)
end


function plot_evolution(gens; stride=1)
    g_len = length(gens[1][1].genome)

    min_fit = minimum(map(x->x[end].fitness, gens))
    max_fit = maximum(map(x->x[1].fitness, gens))

    xl = (0,1)
    yl = (min_fit, max_fit)

    fit_arrays = [map(x->x.fitness, gens[i]) for i in 1:length(gens)]
    dens_arrays = [map(x->sum(x.genome)/g_len, gens[i]) for i in 1:length(gens)]

    p = scatterplot(dens_arrays[1], fit_arrays[1], xlim=xl, ylim=yl, name="gen 1")
    for j in 1:stride:length(gens)
        if j>1
            scatterplot!(p, dens_arrays[j], fit_arrays[j], name="gen $j")
        end
    end
    p
end
