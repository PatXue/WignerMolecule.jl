import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-left"
tm.init_type = :afm_afe

afm_bias(x, y) = [0, 0, (-1)^(x + div(y,2))]
tm.bias = afm_bias
bias_type = typeof(afm_bias)
tm.B = 0.0
tm.init_B = 1.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 7)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.0025:0.0025:0.05
Ls = [24, 48, 96, 120]
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 100000
    tm.thermalization = 50000 * div(L, 24)
    tm.binsize = div(tm.sweeps, 200)
    for T in Ts
        tm.T = T
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)