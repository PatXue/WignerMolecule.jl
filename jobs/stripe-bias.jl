import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using LinearAlgebra
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-bias"

tm.sweeps = 50000
tm.thermalization = 100000
tm.binsize = 500

stripe_bias(B) = (x, _) -> [0, 0, B * (-1)^(div(x, 2))]
tm.bias = stripe_bias(1.0)
bias_type = typeof(tm.bias)

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
tm.init_T = 10
tm.T = 0.1
Ls = [20, 40, 80]
Bs = 0.0:0.05:0.5
for L in Ls
    tm.Lx = tm.Ly = L
    for B in Bs
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        tm.B = B
        tm.bias = stripe_bias(B)
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)