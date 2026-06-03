import Pkg
Pkg.activate(@__DIR__)

include("Expectation.jl")
using .Expectations

using JLD2
using LinearAlgebra
using WignerMolecule

norm2(v) = sum(abs2.(v))

const L = 12
const phases = Dict(
    "stripe" => (45, 5, 20, 6),
    "fm" => (45, 5, 20, 9),
    "afm_fe" => (45, 11, 20, 10),
    "afm_afe" => (45, 11, 20, 7)
)
const spins = Dict(
    "stripe" => (div(L,4)+1, 1),
    "fm" => (1, 1),
    "afm_fe" => (div(L,2)+1, 1),
    "afm_afe" => (div(L,2)+1, div(L,4)+1)
)
const etas = Dict(
    "stripe" => (div(L,2)+1, 1),
    "fm" => (1, 1),
    "afm_fe" => (1, 1),
    "afm_afe" => (1, div(L,2)+1)
)

function sample(name, ord, n; printfreq=100000)
    if name ∉ keys(phases)
        error("Invalid phase name passed: $name")
    end

    raw_params = load_object("all_params.jld2")[phases[name]]
    norm_params = raw_params ./ norm(raw_params)
    wigparams = WignerParams(norm_params...)
    mc = WignerMC{:HighTemp, Nothing}(; Lx=L, Ly=L, wigparams, bias=nothing)

    all_data = load("expectations.jld2")
    avg_energy = [get(all_data, "$name/HH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_scorr = [get(all_data, "$name/sH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηcorr = [get(all_data, "$name/ηH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηzcorr = [get(all_data, "$name/ηzH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_sk = [get(all_data, "$name/skH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηx = [get(all_data, "$name/ηkxH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηy = [get(all_data, "$name/ηkyH$i", Expectation(0,0,0)) for i in 0:ord]
    avg_ηz = [get(all_data, "$name/ηkzH$i", Expectation(0,0,0)) for i in 0:ord]

    for i in 1:n
        randomize!(mc)
        E = total_energy(mc)
        sk = mc.spinks[spins[name]..., :]
        sk_corr = norm2(sk)
        ηk = mc.ηks[etas[name]..., :]
        ηk_corr = norm2(ηk[1:2])    # In-plane η correlation
        ηz_corr = abs2(ηk[3])       # ηz correlation
        energies = [E^i for i in 0:ord]
        avg_energy .= addsample.(avg_energy, E .* energies)
        avg_scorr .= addsample.(avg_scorr, sk_corr .* energies)
        avg_ηcorr .= addsample.(avg_ηcorr, ηk_corr .* energies)
        avg_ηzcorr .= addsample.(avg_ηzcorr, ηz_corr .* energies)
        avg_sk .= addsample.(avg_sk, real(sk[1]) .* energies)
        avg_ηx .= addsample.(avg_ηx, real(ηk[1]) .* energies)
        avg_ηy .= addsample.(avg_ηy, real(ηk[2]) .* energies)
        avg_ηz .= addsample.(avg_ηz, real(ηk[3]) .* energies)

        if i % printfreq == 0
            println("Sample #$i completed")
        end
    end

    for i in 0:ord
        all_data["$name/HH$i"] = avg_energy[i+1]
        all_data["$name/sH$i"] = avg_scorr[i+1]
        all_data["$name/ηH$i"] = avg_ηcorr[i+1]
        all_data["$name/ηzH$i"] = avg_ηzcorr[i+1]
        all_data["$name/skH$i"] = avg_sk[i+1]
        all_data["$name/ηkxH$i"] = avg_ηx[i+1]
        all_data["$name/ηkyH$i"] = avg_ηy[i+1]
        all_data["$name/ηkzH$i"] = avg_ηz[i+1]
        println("Order $i:")
        println("H^$(i+1): $(avg_energy[i+1])")
        println("sH^$i: $(avg_scorr[i+1])")
        println("ηH^$i: $(avg_ηcorr[i+1])")
        println("ηzH^$i: $(avg_ηzcorr[i+1])")
        println("skH$i: $(avg_sk[i+1])")
        println("ηkxH$i: $(avg_ηx[i+1])")
        println("ηkyH$i: $(avg_ηy[i+1])")
        println("ηkzH$i: $(avg_ηz[i+1])")
        println("")
    end
    save("expectations.jld2", all_data)
end

if "-r" ∈ ARGS
    all_data = load("expectations.jld2")
    d = filter(p->!startswith(p.first, ARGS[1]), all_data)
    save("expectations.jld2", d)
end
sample(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]))