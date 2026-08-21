import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "fm-rand"
tm.init_type = :rand
tm.algtype = :Heatbath
tm.init_T = 1.0
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 250

tm.wigparams = WignerParams("all_params.jld2", 5, 9)
Ts = 0.02:0.02:0.2
Ls = [24, 48, 72]
for L in Ls
    tm.Lx = tm.Ly = L
    tm.corr_rad = 2
    for T in Ts
        tm.T = T
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
