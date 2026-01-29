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

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :stripe

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
tm.bias = stripe_bias
bias_type = typeof(stripe_bias)
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
tm.init_B = 10.0
tm.Lx = tm.Ly = 20
Ts = [0.5, 0.6, 0.7]
Bs = 0.0:0.1:1.0
for T in Ts
    tm.T = T
    for B in Bs
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        tm.B = B
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)