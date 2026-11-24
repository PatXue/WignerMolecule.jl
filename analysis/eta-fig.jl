import Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using Carlo.ResultTools
using CarloAnalysis
using DataFrames
using JLD2
using LinearAlgebra
using NaturalNeighbours
using Statistics
set_theme!(theme_latexfonts())

inch = 96
pt = 4/3
cm = inch / 2.54
width = 8.5cm

function linfit(xs, ys; weight=true)
    A = ones(length(xs), 2)
    A[:, 2] = xs
    if weight
        W = Diagonal(getfield.(ys, :err) .^ (-2))
    else
        W = Diagonal(ones(length(ys)))
    end
    return (A' * W * A) \ (A' * W * ys)
end
# Extrapolate to infinite system size
function extrapolate(Ls, ys)
    β = linfit(Ls .^ (-2), ys, weight=true)
    return β[1]
end

checkeq(x) = y -> x .== y

results = JobResult("../eta-jobs", "full-diagram")
vars = [:sk_corr_fm, :chik_corr_Γ, :chi2avg]
# df_extrap = combine(
#     groupby(results.data, [:Jz, :T]),
#     [[:Lx, x] => extrapolate => x for x in vars]...,
# )
df_extrap = results.data
Ts = unique(df_extrap[:,:T])
Jps = hcat([unique(df[:, :Jp]) for df in groupby(df_extrap, :Jz)]...)
Jzs = collect(hcat([unique(df_extrap[:, :Jz]) for _ in 1:size(Jps, 1)]...)')
coords = vec(collect(Iterators.zip(Jps, Jzs)))

bresults = JobResult("../eta-jobs", "full-diagram-broad")
bTs = unique(bresults.data[:,:T])
bJps = unique(bresults.data[:, :Jp])
bJzs = unique(bresults.data[:, :Jz])

fig = Figure(size=(2width, 0.87width), figure_padding=3, fontsize=9pt)

subf = fig[1,1] = GridLayout()
ax = Axis(subf[1,1], title=L"Order Param. Mag. $(T=0.8)$",
    xlabel=L"J_+", ylabel=L"J_z", aspect=DataAspect(),
    limits=((0,2), (0,2)),
    xgridvisible=false, ygridvisible=false,
    xlabelpadding=0, ylabelpadding=0,
    xticksize=2, yticksize=2,
    xticks = [0, 1, 2], yticks = [0, 1, 2],
    backgroundcolor=:snow3
)

phasetonum = Dict(:none => NaN, :stripe => 1, :fm => 2, :fe => 3)
phases = [:stripe, :fm, :fe]
for i in 1:3
    colormap = [Makie.wong_colors(0)[i], Makie.wong_colors(1)[i]]
    vals = broadcast(((Jp, Jz),) -> begin
        data = subset(bresults.data, :Jp => checkeq(Jp), :Jz => checkeq(Jz))
        return data[1, Symbol(:sk_corr_, phases[i])].val
    end, Iterators.product(bJps, bJzs))
    heatmap!(ax, bJps, bJzs, vals; colormap, colorrange=(0,1))
end
linesegments!(
    ax, [(0,1), (1,1), (1,0), (1,1), (1,1), (1.75,2)],
    linestyle=:dash, linewidth=0.5pt, color=:red
)

Legend(subf[2,1],
    [PolyElement(; color) for color in Makie.wong_colors()[1:3]],
    ["Stripe", "FM", "FE"],
    patchsize=(5,5),
    orientation=:horizontal,
    tellwidth=true, tellheight=true,
    valign=:top,
)
rowsize!(subf, 2, 0.1width)
rowgap!(subf, 5)
Label(fig[1,1,TopLeft()], "(a)")

cols = 3
subf = fig[1,2] = GridLayout()
df = df_extrap
# allchis = sqrt.(abs.(getfield.(df[:, :chik_corr_Γ] - df[:, :chi2avg], :val) ./ df[:,:Lx].^2))
# minchi = round(minimum(allchis), sigdigits=1)
# maxchi = round(maximum(allchis), sigdigits=1)
for i in 1:3
    pos = [cld(i, cols), mod1(i, cols)] + [0,1]
    T = Ts[i]
    ax1 = Axis(
        subf[pos...], title=L"T = %$T", aspect=DataAspect(),
        limits=((0,2), (0,2)),
        xgridvisible=false, ygridvisible=false,
        xlabel = L"J_+", ylabel=L"J_z",
        xlabelvisible=false, ylabelvisible=pos[2] == 2,
        xticksvisible=false, xticklabelsvisible=false,
        yticksvisible=pos[2] == 1, yticklabelsvisible=pos[2] == 2,
        xlabelpadding=0, ylabelpadding=0,
        xticks=[0,1,2], yticks=[1,2],
        xticksize=2, yticksize=2,
    )
    ax2 = Axis(
        subf[pos+[1,0]...], aspect=DataAspect(),
        limits=((0,2), (0,2)),
        xgridvisible=false, ygridvisible=false,
        xlabel = L"J_+", ylabel=L"J_z",
        xlabelvisible=true, ylabelvisible=pos[2] == 2,
        yticksvisible=pos[2] == 1, yticklabelsvisible=pos[2] == 2,
        xlabelpadding=0, ylabelpadding=0,
        xticks=[0,1,2], yticks=[1,2],
        xticksize=2, yticksize=2,
    )

    eq = setX -> (x -> x .≈ setX)
    local df = subset(df_extrap, :T => checkeq(Ts[i]))
    chis = broadcast((Jp, Jz) -> begin
        df = subset(df_extrap, :T => checkeq(Ts[i]), :Jz => checkeq(Jz), :Jp => checkeq(Jp))
        interchi = getfield(real(df[1, :chik_corr_Γ] - df[1, :chi2avg]), :val) * df[1,:Lx]^(-2)
        return sqrt(abs(interchi))
    end, Jps, Jzs)
    fms = broadcast((Jp, Jz) -> begin
        df = subset(df_extrap, :T => checkeq(Ts[i]), :Jz => checkeq(Jz), :Jp => checkeq(Jp))
        return getfield(df[1, :sk_corr_fm], :val)
    end, Jps, Jzs)

    hm1 = contourf!(ax1, Jps, Jzs, chis, levels=range(0, 0.0305, 20))
    hm2 = contourf!(ax2, Jps, Jzs, fms, levels=range(0, 1, 20))
    linesegments!(
        ax1, [(0,1), (1,1), (1,0), (1,1), (1,1), (1.75,2)],
        linestyle=:dash, linewidth=0.5pt, color=:red
    )
    linesegments!(
        ax2, [(0,1), (1,1), (1,0), (1,1), (1,1), (1.75,2)],
        linestyle=:dash, linewidth=0.5pt, color=:red
    )
    tightlimits!(ax1)
    tightlimits!(ax2)

    if i == cols
        Colorbar(
            subf[1, cols+2], hm1,
            label=L"\tilde\chi",
            tellheight=false,
            halign=:left,
            ticks=[0, 0.01, 0.02, 0.03],
            ticksize=2, labelpadding=1,
            size = 8, height=0.3width
        )
        Colorbar(
            subf[2, cols+2], hm2,
            label=L"FM Order $\eta^{zz}(\Gamma)$",
            tellheight=false,
            halign=:left,
            ticksize=2, labelpadding=1,
            size = 8, height=0.3width
        )
        rowgap!(subf, 5)
        colgap!(subf, 5)
    end
end

Label(subf[:,1,Top()], "(b)")
colsize!(subf, 1, 0.1pt)

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
colsize!(fig.layout, 1, Relative(0.3))
save("eta-plots/main_fig.pdf", fig)
fig