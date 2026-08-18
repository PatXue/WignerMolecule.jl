import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "stripe"
tm.init_type = :stripe
tm.wigparams = WignerParams("all_params.jld2", 5, 6)

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 250
Ts = 0.02:0.02:0.2
Ls = [48]
fugs = [0.01, 0.05, 0.15, 0.25, 0.5, 1.0]
for (T, L, fug) in Iterators.product(Ts, Ls, fugs)
    tm.fug = fug
    tm.Lx = tm.Ly = L
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)