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
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

tm.wigparams = WignerParams("all_params.jld2", 11, 7)
Ts = range(0.0325, 0.0375, 10)
Ls = [48, 72, 96]
tm.parallel_tempering = (
    mc = WignerMC,
    parameter = :T,
    values = Ts,
    interval = 5
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