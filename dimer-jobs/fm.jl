import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm"
tm.init_type = :fm
tm.wigparams = WignerParams("all_params.jld2", 5, 9)

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 250
Ts = 0.02:0.02:0.2
Ls = [24]
fugs = 0.05:0.05:1.0
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