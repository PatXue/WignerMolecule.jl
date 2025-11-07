import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "stripe"

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :stripe

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 5, 20, 6)]...)
Ts = 0.0:0.25:2.0
for T in Ts
    tm.T = max(T, 0.01)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)