import Pkg
Pkg.activate("..")

using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

fm_bias(_, _) = [0, 0, 1]
bias = fm_bias
bias_type = typeof(fm_bias)

Lx = Ly = 40
B = 0.0

raw_params = load_object("all_params.jld2")[(45, 5, 20, 9)]
norm_params = raw_params ./ norm(raw_params)
wigparams = WignerParams(norm_params...)

mc = WignerMC{:None, bias_type}(; wigparams, bias, Lx, Ly)
for xi in 1:20
    init_skyrm!(mc, xi)
    E = total_energy(mc, B)
    println("Size $xi: $E")
end
