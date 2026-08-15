import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "er7"
tm.algtype = :Heatbath
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :rand
tm.init_T = 5.0

Ls = [48]
Ts = [0.5, 0.75, 1.0, 2.0]
ams = 4:11
for (am, T, L) in Iterators.product(ams, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.corr_rad = div(L, 12)
    tm.am = am
    tm.wigparams = WignerParams("all_params.jld2", 7, am, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)