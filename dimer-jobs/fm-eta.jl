import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm-eta"
tm.algtype = :EtaOnly
tm.wigparams = WignerParams("all_params.jld2", 5, 9)
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.fug = 1

n = 20
Ls = [72]
inits = [:fm, :nematic, :nematic2, :zigzag, :rand]
Ts = [0.005, 0.01]
for (T, init, L) in Iterators.product(Ts, inits, Ls)
    tm.Lx = tm.Ly = L
    tm.T = T
    tm.init_T = 0.1
    tm.init_type = init
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)
