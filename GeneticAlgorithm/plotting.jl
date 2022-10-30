#=
using UnicodePlots
=#

function plot_fit_dens(pop)
    g_len = length(pop[1].genome)
    fitss = [a.fitness for a in pop]
    dens = [sum(a.genome)/g_len for a in pop]

    xl = (0,1)
    yl = (0, 1.15*maximum(fitss))

    scatterplot(dens, fitss, xlim=xl, ylim=yl, xlabel="AV link density", ylabel="Fitness")
end


function plot_evolution(gens; stride=1)
    g_len = length(gens[1][1].genome)

    min_fit = minimum(map(x->x[end].fitness, gens))
    max_fit = maximum(map(x->x[1].fitness, gens))

    xl = (0,1)
    yl = (min_fit, max_fit)

    fit_arrays = [map(x->x.fitness, gens[i]) for i in 1:length(gens)]
    dens_arrays = [map(x->sum(x.genome)/g_len, gens[i]) for i in 1:length(gens)]

    p = scatterplot(dens_arrays[1], fit_arrays[1], xlim=xl, ylim=yl, name="gen 1", xlabel="AV link density", ylabel="Fitness")
    for j in 1:stride:length(gens)
        if j>1
            scatterplot!(p, dens_arrays[j], fit_arrays[j], name="gen $j")
        end
    end
    p
end

"""
Plotting in terminal of network with highlighting of links
for AV exclusive use.

Based on `term_plot` in TrafficNetworks2
"""
function term_link_plot(rn::RoadNetwork, a::T, map_dict::Dict{Int64,Int64}) where {T<:AbstractIndividual}
    node_positions = rn.node_params[:pos]
    edge_num = ne(rn.g)

    genome_ones = findall(x->x==1, a.genome)
    av_edges = [map_dict[i] for i in genome_ones]

    edge_coords = Array{NTuple{4,Float64},1}(undef, edge_num)
    for (i, ed) in enumerate(edges(rn.g))
        x_src, y_src = node_positions[ed.src,:]
        x_dst, y_dst = node_positions[ed.dst,:]
        edge_coords[i] = (x_src, y_src, x_dst, y_dst)
    end

    canvas = BrailleCanvas(60,
                           25,
                           origin_x = 0.0,
                           origin_y = 0.0,
                           width = 1.0,
                           height = 1.0)
    for i in 1:edge_num
        if i ∉ av_edges
            lines!(canvas, edge_coords[i]...)
        end

    end

    for i in 1:edge_num
        if i in av_edges
            lines!(canvas, edge_coords[i]..., color=:red)
        end
    end

    canvas
end
