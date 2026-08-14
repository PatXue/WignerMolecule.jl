import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-vbs"
tm.init_type = :vbs
tm.wigparams = WignerParams("all_params.jld2", 11, 10)

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 500
Ls = [48]
fugs = [0.1, 0.2, 0.5, 1.0, 5.0, 10.0]
Ts = 0.05:0.05:0.5
for (T, fug, L) in Iterators.product(Ts, fugs, Ls)
    tm.Lx = tm.Ly = L
    tm.T = T
    tm.fug = fug
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)