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

ord = 4
all_data = load("expectations.jld2")
avg_energy = [get(all_data, "stripe/H$i", Expectation(0,0,0)) for i in 1:ord]

n = 10^5
for _ in 1:n
    randomize!(mc)
    E = total_energy(mc)
    sk = mc.spinks[div(L, 4)+1, 1, :]
    ηk = mc.ηks[div(Lx, 2)+1, 1, :]
    avg_energy .= addsample.(avg_energy, [E^i for i in 1:ord])
end

for i in 1:ord
    all_data["stripe/H$i"] = avg_energy[i]
    println("H^$i: $(avg_energy[i])")
end
save("expectations.jld2", all_data)
