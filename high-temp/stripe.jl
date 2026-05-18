import Pkg
Pkg.activate(@__DIR__)

include("Expectation.jl")
using .Expectations

using JLD2
using LinearAlgebra
using WignerMolecule

norm2(v) = sum(abs2.(v))

Lx = Ly = L = 8
raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
wigparams = WignerParams(norm_params...)
mc = WignerMC{:HighTemp, Nothing}(; Lx, Ly, wigparams, bias=nothing)

avg_energy = Expectation(0, 0, 0)
n = 10^3
ord = 4
for _ in 1:n
    randomize!(mc)
    E = total_energy(mc)
    sk = mc.spinks[div(L, 4)+1, 1, :]
    ηk = mc.ηks[div(Lx, 2)+1, 1, :]
    global avg_energy = addsample(avg_energy, E)
end

println(avg_energy)
