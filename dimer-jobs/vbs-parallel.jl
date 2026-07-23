import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-parallel"
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

tm.sweeps = 50000
tm.binsize = 250
Ts = 0.01:0.01:0.15
Ls = [24]
tm.parallel_tempering = (
    mc = DimerMC,
    parameter = :T,
    values = Ts,
    interval = 10
)
for L in Ls
    tm.thermalization = 100000
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