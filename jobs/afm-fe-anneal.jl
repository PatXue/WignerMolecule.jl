import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-anneal"

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = div(tm.sweeps, 100)
tm.init_type = :afm_fe

afm_bias(x, _) = [0, 0, (-1)^x]
tm.bias = afm_bias
bias_type = typeof(afm_bias)
tm.B = 0.0
tm.init_B = 10.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.05:0.05:0.6
Ls = [20]
for L in Ls
    tm.Lx = tm.Ly = L
    for T in Ts
        tm.T = max(T, 0.01)
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)