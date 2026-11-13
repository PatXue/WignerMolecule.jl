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
tm.thermalization = 400000
tm.binsize = 1000

n = 20
Ls = [72]
Ts = [0.02, 0.04, 0.06, 0.08]
tm.betamus = range(-30, 0, 20)
for (L, T) in Iterators.product(Ls, Ts)
    tm.Lx = tm.Ly = L
    tm.parallel_tempering = (
        mc = DimerMC,
        parameter = :fug,
        values = exp.(tm.betamus),
        interval = 1
    )
    tm.T = T
    tm.init_T = T + 0.1
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
