import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-anneal"
tm.init_type = :stripe
tm.corr_rad = 3
tm.wigparams = WignerParams("all_params.jld2", 5, 6)

tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000
Ts = [0.01, 0.02, 0.04, 0.08, 0.1]
Ls = [60]
βμs = range(-7, 2, 20)
for (L, βμ, T) in Iterators.product(Ls, βμs, Ts)
    tm.Lx = tm.Ly = L
    tm.betamu = βμ
    tm.fug = exp(-βμ)
    tm.T = T
    tm.init_T = T + 0.2
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)