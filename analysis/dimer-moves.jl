import Pkg
Pkg.activate(@__DIR__)

using CairoMakie
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

fig = Figure(size=(width, width), figure_padding=1)
markersize = pt
posns = [[0,0], [1,0], [0,1], [1,1]]
plaquette = idxtoreal.(posns)
hollow_circle = BezierPath([
    MoveTo(Point(1, 0)),
    EllipticalArc(Point(0, 0), 1, 1, 0, 0, 2pi),
    MoveTo(Point(0.9, 0.9)),
    EllipticalArc(Point(0, 0), 0.9, 0.9, 0, 0, 2pi),
    ClosePath(),
])

for pos in posns
    pos += [1,1]
    subf = fig[pos...] = GridLayout()
    ax = Axis(subf[1,1], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    scatter!(ax, plaquette; marker=hollow_circle, markersize, color=:black)
    ax = Axis(subf[3,1], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    scatter!(ax, plaquette; marker=hollow_circle, markersize, color=:black)
end

save("dimer-plots/dimer_moves.pdf", fig)