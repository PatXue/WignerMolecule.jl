import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools
using JLD2

tm = TaskMaker()
jobname = "stripe-eta"

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 0
tm.binsize = 500
tm.init_type = :stripe_s

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.0:0.5:5.0
for T in Ts
    # tm.thermalization = T < 0.3 ? 40000 : 20000
    tm.T = max(T, 0.01)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo(jobname, WignerMC{:Metropolis_η, Nothing},
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)