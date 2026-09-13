import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-crit"
tm.init_type = :stripe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 2000

tm.wigparams = WignerParams("all_params.jld2", 5, 6)
Ls = [48, 60, 72]
Ts = range(0.08, 0.0812, 12)
tm.Ts = Ts
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