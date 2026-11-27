import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-anneal"
tm.init_type = :stripe
tm.corr_rad = 2
tm.wigparams = WignerParams("all_params.jld2", 5, 6)

n = 20
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000
Ls = [72]
Ts = [0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.0075, 0.005, 0.0025]
tm.betamus = range(-20, 20, n)
for (L, T) in Iterators.product(Ls, Ts)
    tm.Lx = tm.Ly = L
    tm.parallel_tempering = (
        mc = DimerMC,
        parameter = :fug,
        values = exp.(tm.betamus),
        interval = 1
    )
    tm.T = T
    tm.init_T = 0.08
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
