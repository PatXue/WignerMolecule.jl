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

ord = 5
all_data = load("expectations.jld2")
avg_energy = [get(all_data, "stripe/HH$i", Expectation(0,0,0)) for i in 0:ord]
avg_sk = [get(all_data, "stripe/sH$i", Expectation(0,0,0)) for i in 0:ord]
avg_ηk = [get(all_data, "stripe/ηH$i", Expectation(0,0,0)) for i in 0:ord]
avg_ηz = [get(all_data, "stripe/ηzH$i", Expectation(0,0,0)) for i in 0:ord]

n = 10^6
for _ in 1:n
    randomize!(mc)
    E = total_energy(mc)
    sk = mc.spinks[div(L, 4)+1, 1, :]
    sk_corr = norm2(sk)
    ηk = mc.ηks[div(Lx, 2)+1, 1, :]
    ηk_corr = norm2(ηk[1:2])    # In-plane η correlation
    ηz_corr = abs2(ηk[3])       # ηz correlation
    energies = [E^i for i in 0:ord]
    avg_energy .= addsample.(avg_energy, E .* energies)
    avg_sk .= addsample.(avg_sk, sk_corr .* energies)
    avg_ηk .= addsample.(avg_ηk, ηk_corr .* energies)
    avg_ηz .= addsample.(avg_ηz, ηz_corr .* energies)
end

for i in 0:ord
    all_data["stripe/HH$i"] = avg_energy[i+1]
    all_data["stripe/sH$i"] = avg_sk[i+1]
    all_data["stripe/ηH$i"] = avg_ηk[i+1]
    all_data["stripe/ηzH$i"] = avg_ηz[i+1]
    println("H^$(i+1): $(avg_energy[i+1])")
    println("sH^$i: $(avg_sk[i+1])")
    println("ηH^$i: $(avg_ηk[i+1])")
    println("ηzH^$i: $(avg_ηz[i+1])")
end
save("expectations.jld2", all_data)
