import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-hight"
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

tm.sweeps = 50000
tm.binsize = 250
Ts = 0.2:0.2:4.0
Ls = [24]
for (T, L) in Iterators.product(Ts, Ls)
    tm.thermalization = 50000
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