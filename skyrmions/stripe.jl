import Pkg
Pkg.activate("..")

using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
bias = stripe_bias
bias_type = typeof(stripe_bias)

Lx = Ly = 40
B = 0.0

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
wigparams = WignerParams(norm_params...)

mc = WignerMC{:None, bias_type}(; wigparams, bias, Lx, Ly)
for xi in 1:20
    init_skyrm!(mc, xi, etatype=:stripe)
    for I in eachindex(mc.spins)
        x, y = Tuple(I)
        mc.spins[I] = mc.spins[I] .* (-1)^(div(x,2))
    end
    E = total_energy(mc, B)
    println("Size $(lpad(string(xi), 2, '0')): $E")
end
