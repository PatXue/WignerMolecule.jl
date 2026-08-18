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
tm.algtype = :Heatbath

tm.wigparams = WignerParams("all_params.jld2", 5, 6)
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
Ls = [24, 48, 72]
Ts = 0.01:0.01:0.15
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