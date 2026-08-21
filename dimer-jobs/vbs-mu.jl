import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-mu"
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

tm.sweeps = 50000
tm.thermalization = 100000
tm.binsize = 250
Ts = 0.02:0.02:0.2
Ls = [48]
mus = -0.4:0.2:0.4
for (T, L, mu) in Iterators.product(Ts, Ls, mus)
    tm.Lx = tm.Ly = L
    tm.T = T
    tm.mu = mu
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)