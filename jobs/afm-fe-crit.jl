import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-crit"
tm.estimate_covariance = true
tm.init_type = :afm_fe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 400000
tm.binsize = 1000

n = 20
tm.wigparams = WignerParams("all_params.jld2", 11, 10)
Ls = [12, 16, 20, 36, 48]
for L in Ls
    global Ts = range(0.0512, 0.052 + (0.004 * (12/L)^2), n)
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
