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
tm.estimate_covariance = true
tm.meastype = :ScalarOnly
tm.init_type = :stripe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 400000
tm.binsize = 10

n = 20
tm.wigparams = WignerParams("all_params.jld2", 5, 6)
Ls = [48, 60, 72, 84, 96]
Tmin = 0.08025
for L in Ls
    Ts = range(Tmin, Tmin + (0.0806 - Tmin) * (48/L)^2.1, n)
    tm.init_T = 0.0806
    tm.parallel_tempering = (
        mc = WignerMC,
        parameter = :T,
        values = Ts,
        interval = 1
    )
    tm.Lx = tm.Ly = L
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
