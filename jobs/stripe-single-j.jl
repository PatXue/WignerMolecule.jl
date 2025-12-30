import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-single-j"

L = 40
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :stripe
Ts = 0.0:0.2:2.0
min_T = 0.05

tm.J = :EAM_Re
tm.wigparams = WignerParams(0, 0, 0, 1.0, 0, 0, 0, 0)
for T in Ts
    tm.T = max(T, min_T)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end
tm.J = :EAM_Im
tm.wigparams = WignerParams(0, 0, 0, 1.0im, 0, 0, 0, 0)
for T in Ts
    tm.T = max(T, min_T)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

tm.J = :PM_Re
tm.wigparams = WignerParams(0, 0, 0, 0, 1.0, 0, 0, 0)
for T in Ts
    tm.T = max(T, min_T)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end
tm.J = :PM_Im
tm.wigparams = WignerParams(0, 0, 0, 0, 1.0im, 0, 0, 0)
for T in Ts
    tm.T = max(T, min_T)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis_η, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)