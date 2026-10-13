import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-anneal"
tm.init_type = :stripe
tm.algtype = :ClusterEta

tm.wigparams = WignerParams("all_params.jld2", 5, 6)
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
Ls = [48, 72, 96]
Ts = range(0.01, 0.1, 10)
for L in Ls
    tm.Lx = tm.Ly = L
    for T in Ts
        tm.T = T
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)