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
tm.init_type = :afm_fe

afm_bias(x, _) = [0, 0, (-1)^x]
tm.bias = afm_bias
bias_type = typeof(afm_bias)
tm.B = 0.0
tm.init_B = 5.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.004:0.004:0.08
Ls = [24, 48, 96]
for (T, L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 24)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)