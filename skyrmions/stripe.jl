import Pkg
Pkg.activate("..")

using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
bias = stripe_bias
bias_type = typeof(stripe_bias)
JSON.lower(f::bias_type) = f(1, 1)

Lx = Ly = 40

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
wigparams = WignerParams(norm_params...)

mc = WignerMC{:None, bias_type}(; wigparams, bias, Lx, Ly)
