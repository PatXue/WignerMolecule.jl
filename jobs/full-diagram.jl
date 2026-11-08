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
tm.corr_rad = 3

tm.init_T = 0.12
Ls = [84]
Ts = [0.007, 0.01, 0.02, 0.04, 0.08, 0.12]
ams = 4:0.5:11
ers = 5:0.5:11
for (am, er, T, L) in Iterators.product(ams, ers, Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.am = am
    tm.er = er
    tm.wigparams = WignerParams("params_interpolated_quarter.jld2", (am, er))
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
