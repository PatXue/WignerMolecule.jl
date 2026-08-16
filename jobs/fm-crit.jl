import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "fm-crit"

tm.algtype = :Heatbath
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
tm.init_type = :const

tm.wigparams = WignerParams("all_params.jld2", 5, 9)
Ts = 0.09:0.0025:0.1175
Ls = [24, 48, 96]
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