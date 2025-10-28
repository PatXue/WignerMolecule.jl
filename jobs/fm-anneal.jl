import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm-anneal"

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 20000
tm.thermalization = 20000
tm.binsize = 100

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 5, 20, 9)]...)
tm.init_T = 10
Ts = 0.5:0.5:7
Ls = [20, 30, 40]
for L in Ls
    tm.Lx = tm.Ly = L
    for T in Ts
        tm.T = max(T, 0.01)
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)