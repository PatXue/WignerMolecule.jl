import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm-vbs"
tm.init_type = :vbs
tm.wigparams = WignerParams("all_params.jld2", 5, 9)

tm.sweeps = 50000
tm.thermalization = 100000
tm.binsize = 250
Ts = 0.01:0.01:0.2
Ls = [24, 48]
for (T, L) in Iterators.product(Ts, Ls)
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