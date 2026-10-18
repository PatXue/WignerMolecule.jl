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
tm.init_type = :afm_afe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 400000
tm.binsize = 1000

tm.init_T = 0.0365
tm.wigparams = WignerParams("all_params.jld2", 11, 7)
Ls = [12, 18, 24]
for L in Ls
    tm.Lx = tm.Ly = L
    Ts = range(0.036, 0.036 + 0.001 * (12/L)^2, 12)
    tm.parallel_tempering = (
        mc = WignerMC,
        parameter = :T,
        values = Ts,
        interval = 1
    )
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = length(Ts)
)
start(job, ARGS)