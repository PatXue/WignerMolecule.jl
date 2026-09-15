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
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
tm.fug = 1.0

Ls = [48]
Ts = [1.0, 2.0, 4.0]
ams = 4:11
ers = 5:11
for (am, er, T, L) in Iterators.product(ams, ers, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.am = am
    tm.er = er
    tm.wigparams = WignerParams("all_params.jld2", er, am, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
