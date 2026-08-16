import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm-anneal"
tm.init_T = 2.0
tm.init_type = :rand
tm.wigparams = WignerParams("all_params.jld2", 5, 9)

tm.sweeps = 50000
tm.thermalization = 100000
tm.binsize = 500
Ts = Iterators.flatten((0.02:0.01:0.09, 0.1:0.025:0.2))
Ls = [48]
fugs = [0.05, 0.5, 1.0, 2.0]
for (T, L, fug) in Iterators.product(Ts, Ls, fugs)
    tm.Lx = tm.Ly = L
    tm.corr_rad = div(L, 24)
    tm.fug = fug
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)