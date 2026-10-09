import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-crit-rand"
tm.init_type = :rand
tm.algtype = :Heatbath
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000

tm.wigparams = WignerParams("all_params.jld2", 5, 6)
Ls = [60, 72, 84]
Ts = collect(range(0.07, 0.08, 12))
tm.Ts = Ts
tm.init_T = 1.0
tm.corr_rad = 2
tm.parallel_tempering = (
    mc = WignerMC,
    parameter = :T,
    values = Ts,
    interval = 1
)
for L in Ls
    tm.Lx = tm.Ly = L
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = length(Ts)
)
start(job, ARGS)