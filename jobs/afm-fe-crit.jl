import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-crit"
tm.init_type = :afm_fe
tm.algtype = :Cluster
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

tm.wigparams = WignerParams("all_params.jld2", 11, 10)
Ts = collect(Iterators.flatten((0.0513:0.0001:0.0523,)))
Ls = [36, 48, 60, 72]
tm.parallel_tempering = (
    mc = WignerMC,
    parameter = :T,
    values = Ts,
    interval = 10
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