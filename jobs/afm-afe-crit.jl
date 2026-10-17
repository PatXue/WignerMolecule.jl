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
tm.sweeps = 400000
tm.thermalization = 400000
tm.binsize = 1000

tm.init_T = 0.0365
tm.wigparams = WignerParams("all_params.jld2", 11, 7)
Ts = collect(range(0.0359, 0.0365, 12))
Ls = [24, 36, 48]
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