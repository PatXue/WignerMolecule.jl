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

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 500
Ts = [0.001, 0.01, 0.02, 0.04, 0.08, 0.12]
Ls = [60]
βμs = range(-4, 4, 20)
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