import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs"
tm.Q = 0.5
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

Ts = 0.02:0.02:0.2
Ls = [24]
for (T, L) in Iterators.product(Ts, Ls)
    tm.sweeps = 50000 * div(L, 24)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
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