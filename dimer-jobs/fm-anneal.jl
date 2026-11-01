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
tm.corr_rad = 2
tm.wigparams = WignerParams("all_params.jld2", 5, 9)

tm.sweeps = 50000
tm.thermalization = 100000
tm.binsize = 500
Ts = range(0.02, 0.15, 5)
Ls = [60]
βμs = range(-1, 3, 10)
for (L, βμ, T) in Iterators.product(Ls, βμs, Ts)
    tm.Lx = tm.Ly = L
    tm.betamu = βμ
    tm.fug = exp(-βμ)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)