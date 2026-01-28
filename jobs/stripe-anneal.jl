import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-anneal"

L = 20
tm.Lx = tm.Ly = L
tm.init_type = :stripe

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
tm.bias = stripe_bias
bias_type = typeof(stripe_bias)
tm.B = 0.0
tm.init_B = 10.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ls = [20, 40, 80]
Ts = 0.1:0.1:1.0
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 100000 * div(L, 20)
    tm.thermalization = 100000 * div(L, 20)
    tm.binsize = 500
    for T in Ts
        tm.T = T
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)