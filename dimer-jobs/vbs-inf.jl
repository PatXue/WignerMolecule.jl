import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-inf"
tm.init_type = :rand
tm.algtype = :Worm
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

tm.sweeps = 10000
tm.thermalization = 10000
tm.binsize = 100
tm.T = Inf
betamus = Iterators.flatten((range(-10, 10, 20), 0))
Ls = [48]
for (betamu, L) in Iterators.product(betamus, Ls)
    tm.Lx = tm.Ly = L
    tm.betamu = betamu
    tm.fug = exp(betamu)
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)