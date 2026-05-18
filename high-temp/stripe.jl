import Pkg
Pkg.activate(@__DIR__)

using JLD2
using WignerMolecule

Lx = Ly = 6
raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
wigparams = WignerParams(norm_params...)
mc = WignerMC{:HighTemp, Nothing}(; Lx, Ly, wigparams, bias=nothing)
