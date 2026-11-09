import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-anneal"
tm.init_T = 2.0
tm.init_type = :rand
tm.algtype = :Worm
tm.corr_rad = 2
tm.wigparams = WignerParams("all_params.jld2", 10, 7)

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 250
Ts = [0.02, 0.05]
Ls = [72]
βμs = range(-10, 10, 20)
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
