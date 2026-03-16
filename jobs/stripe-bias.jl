import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using LinearAlgebra
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-bias"
tm.init_type = :stripe

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
tm.bias = stripe_bias
bias_type = typeof(stripe_bias)
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
tm.init_B = 5.0
Ts = [0.5, 0.6, 0.7]
Bs = 0.0:0.05:0.5
Ls = [20, 40, 80]
for (B, T, L) in Iterators.product(Bs, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 20)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)

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