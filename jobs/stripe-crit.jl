import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-crit"
tm.init_type = :stripe
tm.algtype = :Heatbath

tm.wigparams = WignerParams("all_params.jld2", 5, 6)
tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 1000
Ls = [48, 72, 96]
Ts = 0.08:0.0002:0.082
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