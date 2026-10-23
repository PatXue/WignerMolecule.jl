import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "full-diagram"
tm.init_type = :vbs
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000

Ls = [60]
fugs = [0.1, 0.5, 1.0]
Ts = [0.5, 1.0, 2.0, 4.0]
ams = 4:11
ers = 5:11
for (am, er, T, fug, L) in Iterators.product(ams, ers, Ts, fugs, Ls)
    tm.Lx = tm.Ly = L
    tm.fug = fug
    tm.am = am
    tm.er = er
    tm.wigparams = WignerParams("all_params.jld2", er, am, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo(jobname, DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
