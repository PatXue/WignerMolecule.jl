import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "full-diagram"
tm.algtype = :Heatbath
tm.init_type = :rand
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
tm.corr_rad = 2

tm.init_T = 0.1
Ls = [48]
Ts = [0.005, 0.007, 0.01, 0.02, 0.04, 0.08]
ams = 4:11
ers = 5:11
for (am, er, T, L) in Iterators.product(ams, ers, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.am = am
    tm.er = er
    tm.wigparams = WignerParams("all_params.jld2", er, am)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
