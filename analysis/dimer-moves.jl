import Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using CairoMakie.GeometryBasics
using Carlo.ResultTools
using CarloAnalysis
using DataFrames
using HDF5
using JLD2
using LinearAlgebra
using Statistics
set_theme!(theme_latexfonts())

inch = 96
pt = 4/3
cm = inch / 2.54
width = 8.5cm

idxtoreal(v) = v[1] * [1, 0] + v[2] * [1/2, √3/2]

fig = Figure(size=(width, 0.5width), figure_padding=1, fontsize=9pt)
markersize = 3pt
bondwidth = 2pt
posns = [[0,0], [1,0], [0,1], [1,1]]
labels = ["(a)", "(b)", "(c)", "(d)"]
plaquette = Tuple.(5 .* idxtoreal.(posns))
# hollow_circle = BezierPath([
#     MoveTo(Point(1, 0)),
#     EllipticalArc(Point(0, 0), 1, 1, 0, 0, 2pi),
#     MoveTo(Point(0.9, 0)),
#     EllipticalArc(Point(0, 0), 0.9, 0.9, 0, 0, 2pi),
#     ClosePath(),
# ])

initbonds = [
    [plaquette[2], plaquette[3]],
    [],
    [plaquette[1], plaquette[3]],
    [plaquette[1], plaquette[2], plaquette[3], plaquette[4]],
]
finalbonds = [
    [],
    [plaquette[2], plaquette[3]],
    [plaquette[1], plaquette[2]],
    [plaquette[1], plaquette[3], plaquette[2], plaquette[4]],
]

subfs = []
for i in eachindex(posns, labels, initbonds, finalbonds)
    pos = posns[i] + [1,1]
    subf = fig[pos...] = GridLayout()

    ax = Axis(subf[1,1], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    if length(initbonds[i]) != 0
        linesegments!(ax, initbonds[i]; linewidth=bondwidth, color=Makie.wong_colors()[1])
    end
    lines!(ax, plaquette[2:3]; linewidth=pt, linestyle=:dash, color=:red)
    scatter!(ax, plaquette; marker=:circle, markersize, color=:black)

    ax = Axis(subf[1,3], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    if length(finalbonds[i]) != 0
        linesegments!(ax, finalbonds[i]; linewidth=bondwidth, color=Makie.wong_colors()[1])
    end
    scatter!(ax, plaquette; marker=:circle, markersize, color=:black)

    ax = Axis(subf[1,2], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    arrows2d!(
        ax, [(0,0)], [(1,0)];
        shaftwidth=pt,
        tipwidth=3pt, tiplength=0.1cm
    )

    Label(subf[1,1,TopLeft()], labels[i])

    colsize!(subf, 2, Auto(0.25))
    colgap!(subf, 5)
    rowgap!(subf, 5)
    push!(subfs, subf)
end

colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)
save("dimer-plots/dimer_moves.pdf", fig)