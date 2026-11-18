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

checkeq(x) = y -> x .== y

results = JobResult("../eta-jobs", "full-diagram")
Ts = unique(results.data[:,:T])
Jps = hcat([unique(df[:, :Jp]) for df in groupby(results.data, :Jz)]...)
Jzs = collect(hcat([unique(results.data[:, :Jz]) for _ in 1:size(Jps, 1)]...)')
coords = vec(collect(Iterators.zip(Jps, Jzs)))

bresults = JobResult("../eta-jobs", "full-diagram-broad")
bTs = unique(bresults.data[:,:T])
bJps = unique(bresults.data[:, :Jp])
bJzs = unique(bresults.data[:, :Jz])

fig = Figure(size=(2width, 0.87width), figure_padding=3, fontsize=9pt)

subf = fig[1,1] = GridLayout()
ax = Axis(subf[1,1], title=L"\textbf{(a)} $T=0.8$",
    xlabel=L"J_+", ylabel=L"J_z", aspect=DataAspect(),
    limits=((0,2), (0,2)),
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

cols = 3
subf = fig[1,2] = GridLayout()
df = results.data
allchis = getfield.(real.(df[:, :chik_corr_Γ] - df[:, :chi2avg]) ./ df[:,:Lx].^2, :val)
minchi = round(minimum(allchis), sigdigits=1)
maxchi = round(maximum(allchis), sigdigits=1)
for i in 1:3
    pos = [cld(i, cols), mod1(i, cols)]
    T = Ts[i]
    ax1 = Axis(
        subf[pos...], title=L"\textbf{(%$([:b,:c,:d][i]))} $T = %$T$", aspect=DataAspect(),
        limits=((0,2), (0,2)),
        xlabel = L"J_+", ylabel=L"J_z",
        xlabelvisible=false, ylabelvisible=pos[2] == 1,
        xticksvisible=false, xticklabelsvisible=false,
        yticksvisible=pos[2] == 1, yticklabelsvisible=pos[2] == 1,
        xlabelpadding=0, ylabelpadding=0,
        xticks=[0,1,2], yticks=[1,2],
        xticksize=2, yticksize=2,
    )
    ax2 = Axis(
        subf[pos+[1,0]...], aspect=DataAspect(),
        limits=((0,2), (0,2)),
        xlabel = L"J_+", ylabel=L"J_z",
        xlabelvisible=true, ylabelvisible=pos[2] == 1,
        yticksvisible=pos[2] == 1, yticklabelsvisible=pos[2] == 1,
        xlabelpadding=0, ylabelpadding=0,
        xticks=[0,1,2], yticks=[1,2],
        xticksize=2, yticksize=2,
    )

    eq = setX -> (x -> x .≈ setX)
    local df = subset(results.data, :T => checkeq(Ts[i]))
    chis = broadcast((Jp, Jz) -> begin
        df = subset(results.data, :T => checkeq(Ts[i]), :Jz => checkeq(Jz), :Jp => checkeq(Jp))
        return getfield(real(df[1, :chik_corr_Γ] - df[1, :chi2avg]), :val) * df[1,:Lx]^(-2)
    end, Jps, Jzs)
    fms = broadcast((Jp, Jz) -> begin
        df = subset(results.data, :T => checkeq(Ts[i]), :Jz => checkeq(Jz), :Jp => checkeq(Jp))
        return getfield(df[1, :sk_corr_fm], :val)
    end, Jps, Jzs)

    hm1 = contourf!(ax1, Jps, Jzs, chis, levels=range(minchi, 0.001, 20))
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
            subf[1, cols+1], hm1,
            label=L"\left\langle\chi_\Gamma^2\right\rangle - \left\langle\chi_0^2\right\rangle",
            tellheight=false,
            halign=:left,
            ticks=[0, 0.0005, 0.001],
            ticksize=2, labelpadding=1,
            size = 8, height=0.3width
        )
        Colorbar(
            subf[2, cols+1], hm2,
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

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
save("eta-plots/main_fig.pdf", fig)
fig