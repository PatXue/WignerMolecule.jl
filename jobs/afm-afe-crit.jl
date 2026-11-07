import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-crit"
tm.estimate_covariance = true
tm.meastype = :ScalarOnly
tm.init_type = :afm_afe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 10

n = 20
tm.wigparams = WignerParams("all_params.jld2", 11, 7)
Ls = [12, 18, 24, 30, 36, 48]
for L in Ls
    tm.Lx = tm.Ly = L
    Ts = range(0.0359, 0.0359 + 0.0011 * (24/L)^2, n)
    tm.parallel_tempering = (
        mc = WignerMC,
        parameter = :T,
        values = Ts,
        interval = 1
    )
    task(tm)
end

job = JobInfo(jobname, ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
