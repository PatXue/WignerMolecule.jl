import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using LinearAlgebra
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-bias"

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
tm.init_type = :afm_fe

afm_bias(x, _) = [0, 0, (-1)^x]
tm.bias = afm_bias
bias_type = typeof(afm_bias)
tm.init_B = 10.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
tm.init_T = 10
Ls = [20]
Ts = [0.2, 0.3, 0.4]
Bs = 0.0:0.05:0.5
for (B, T, L) in Iterators.product(Bs, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.T = T
    tm.B = B
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)