import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-anneal"

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 100

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 11, 20, 10)]...)
tm.init_T = 10
Ts = 0.0:0.025:0.3
for T in Ts
    tm.T = max(T, 0.01)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)