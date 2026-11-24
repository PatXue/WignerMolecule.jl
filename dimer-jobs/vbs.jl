import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs"
tm.init_type = :vbs
tm.corr_rad = 2

tm.wigparams = WignerParams("all_params.jld2", 10, 7)
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000

n = 20
Ls = [72]
Ts = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
for (T, L) in Iterators.product(Ts, Ls)
    tm.betamus = range(-50, 10, n)
    tm.Lx = tm.Ly = L
    tm.parallel_tempering = (
        mc = DimerMC,
        parameter = :fug,
        values = exp.(tm.betamus),
        interval = 1
    )
    tm.T = T
    tm.init_T = 0.04
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
