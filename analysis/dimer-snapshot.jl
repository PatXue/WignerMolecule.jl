import Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using Carlo.ResultTools
using CarloAnalysis
using DataFrames
using HDF5
using JLD2
using LinearAlgebra
using NaturalNeighbours
using Statistics
set_theme!(theme_latexfonts())

inch = 96
pt = 4/3
cm = inch / 2.54
width = 8.5cm

checkeq(x) = y -> x .== y
idxtoreal(v) = v[1] * [1, 0] + v[2] * [1/2, √3/2]

function generate_spins!(ax, jobname, task_no)
    task_str = lpad(task_no, 4, "0")
    h5open("../dimer-jobs/$jobname.data/task$task_str/run0001.dump.h5") do file
        rawspins = map(t -> [t[:data][1], t[:data][2]], read(file, "simulation/spins"))
        mod_equiv(pos, posj) = all((pos .- posj) .% size(rawspins) .== (0,0))

        # Bounding parallelogram
        Lx, Ly = size(rawspins)
        poly!(ax, Point2f.(idxtoreal.([(0.5,0.5), (0.5,Ly+0.5), (Lx+0.5,Ly+0.5), (Lx+0.5,0.5)])), alpha=0.75)

        positions::Vector{Vector{Float64}} = []
        monomers::Vector{Point2f} = []
        for I in eachindex(IndexCartesian(), rawspins)
            x,y = Tuple(I)
            pos = [x, y]
            disp = rawspins[I] - pos
            for a in ((1,0), (-1,1), (0,-1), (-1,0), (1,-1), (0,1))
                if mod_equiv(a, disp)
                    push!(positions, idxtoreal(pos))
                    push!(positions, idxtoreal(pos .+ a))
                end
            end
            if mod_equiv((0,0), disp)
                push!(monomers, Point2f(idxtoreal(pos)...))
            end
        end
        linesegments!(
            ax, getindex.(positions, 1), getindex.(positions, 2),
            color="purple", linewidth=0.5pt
        )
        scatter!(ax, monomers, color="orange", markersize=pt)
    end
end

function generate_spins!(ax, jobname, task_no, mc_no)
    task_str = lpad(task_no, 4, "0")
    h5open("../dimer-jobs/$jobname.data/task$task_str/run0001.dump.h5") do file
        mc_no = findfirst(x -> x == mc_no, read(file["simulation/chain_permutation"]))
        rawspins = map(
            t -> [t[:data][1], t[:data][2]],
            read(file, "simulation/child_mcs/$mc_no/spins")
        )
        mod_equiv(pos, posj) = all((pos .- posj) .% size(rawspins) .== (0,0))

        # Bounding parallelogram
        Lx, Ly = size(rawspins)
        poly!(ax, Point2f.(idxtoreal.([(0.5,0.5), (0.5,Ly+0.5), (Lx+0.5,Ly+0.5), (Lx+0.5,0.5)])), alpha=0.75)

        positions::Vector{Vector{Float64}} = []
        monomers::Vector{Point2f} = []
        for I in eachindex(IndexCartesian(), rawspins)
            x,y = Tuple(I)
            pos = [x, y]
            disp = rawspins[I] - pos
            for a in ((1,0), (-1,1), (0,-1), (-1,0), (1,-1), (0,1))
                if mod_equiv(a, disp)
                    push!(positions, idxtoreal(pos))
                    push!(positions, idxtoreal(pos .+ a))
                end
            end
            if mod_equiv((0,0), disp)
                push!(monomers, Point2f(idxtoreal(pos)...))
            end
        end
        linesegments!(
            ax, getindex.(positions, 1), getindex.(positions, 2),
            color="purple", linewidth=0.5pt
        )
        scatter!(ax, monomers, color="orange", markersize=pt)
    end
end

jobname = ARGS[1]
j = parse(Int, ARGS[2])
try
    mkdir("dimer-plots/$(jobname)-plots")
catch
end
for i in 1:length(readdir("../dimer-jobs/$jobname.data"))
    fig = Figure(size=(0.5width, 0.3width), figure_padding=1)
    ax = Axis(fig[1,1], aspect=DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)
    generate_spins!(ax, jobname, i, j)
    save("dimer-plots/$(jobname)-plots/$(jobname)$i-$j.pdf", fig)
end