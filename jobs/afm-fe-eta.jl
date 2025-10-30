import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools
using JLD2

tm = TaskMaker()
jobname = "afm-fe-eta"

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 30000
tm.thermalization = 0
tm.binsize = 500
tm.init_type = :afm_fe_s

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 11, 20, 10)]...)
Ts = 0.0:1.0:20.0
for T in Ts
    # tm.thermalization = T < 0.3 ? 40000 : 20000
    tm.T = max(T, 0.01)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo(jobname, WignerMC{:Metropolis_η},
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)