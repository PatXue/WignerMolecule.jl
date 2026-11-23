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

fig = Figure(size=(width, 0.5width), figure_padding=1)
markersize = 2pt
posns = [[0,0], [1,0], [0,1], [1,1]]
plaquette = Tuple.(idxtoreal.(posns))
# hollow_circle = BezierPath([
#     MoveTo(Point(1, 0)),
#     EllipticalArc(Point(0, 0), 1, 1, 0, 0, 2pi),
#     MoveTo(Point(0.9, 0)),
#     EllipticalArc(Point(0, 0), 0.9, 0.9, 0, 0, 2pi),
#     ClosePath(),
# ])

subfs = []
for i in 1:length(posns)
    pos = posns[i] + [1,1]
    subf = fig[pos...] = GridLayout()

    ax = Axis(subf[1,1], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    scatter!(ax, plaquette; marker=:circle, markersize, color=:black)
    ax = Axis(subf[1,2], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    arrows2d!(
        ax, [(0,0)], [(1,0)];
        shaftwidth=0.5pt,
        tipwidth=2pt, tiplength=0.1cm
    )
    ax = Axis(subf[1,3], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    scatter!(ax, plaquette; marker=:circle, markersize, color=:black)

    colsize!(subf, 2, Auto(0.5))
    colgap!(subf, 5)
    rowgap!(subf, 5)
    push!(subfs, subf)
end

save("dimer-plots/dimer_moves.pdf", fig)