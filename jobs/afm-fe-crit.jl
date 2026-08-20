import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-crit"
tm.init_type = :afm_fe
tm.algtype = :Heatbath

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

tm.wigparams = WignerParams("all_params.jld2", 11, 10)
Ts = 0.008:0.008:0.08
Ls = [24, 48, 72]
for (T, L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)