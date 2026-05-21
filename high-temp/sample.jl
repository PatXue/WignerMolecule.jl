import Pkg
Pkg.activate(@__DIR__)

include("Expectation.jl")
using .Expectations

using JLD2
using LinearAlgebra
using WignerMolecule

norm2(v) = sum(abs2.(v))

const phases = Dict(
    "stripe" => (45, 5, 20, 6),
    "fm" => (45, 5, 20, 9),
    "afm_fe" => (45, 11, 20, 10),
    "afm_afe" => (45, 11, 20, 7)
)

function sample(name, ord, n)
    if name ∉ keys(phases)
        error("Invalid phase name passed: $name")
    end

    Lx = Ly = L = 8
    raw_params = load_object("all_params.jld2")[phases[name]]
    norm_params = raw_params ./ norm(raw_params)
    wigparams = WignerParams(norm_params...)
    mc = WignerMC{:HighTemp, Nothing}(; Lx, Ly, wigparams, bias=nothing)

    all_data = load("expectations.jld2")
    avg_energy = [get(all_data, "$name/HH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_sk = [get(all_data, "$name/sH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηk = [get(all_data, "$name/ηH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηz = [get(all_data, "$name/ηzH$i", Expectation(0,0,0)) for i in 0:ord]

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
        all_data["$name/HH$i"] = avg_energy[i+1]
        all_data["$name/sH$i"] = avg_sk[i+1]
        all_data["$name/ηH$i"] = avg_ηk[i+1]
        all_data["$name/ηzH$i"] = avg_ηz[i+1]
        println("H^$(i+1): $(avg_energy[i+1])")
        println("sH^$i: $(avg_sk[i+1])")
        println("ηH^$i: $(avg_ηk[i+1])")
        println("ηzH^$i: $(avg_ηz[i+1])")
    end
    save("expectations.jld2", all_data)
end

sample(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]))