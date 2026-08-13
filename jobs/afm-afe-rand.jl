import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-rand"
tm.algtype = :Heatbath
tm.init_type = :rand
tm.init_T = 0.5

tm.wigparams = WignerParams("all_params.jld2", 11, 7)
Ts = 0.0025:0.0025:0.05
Ls = [24, 48, 96]
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
for (T, L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.corr_rad = div(L, 12)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)