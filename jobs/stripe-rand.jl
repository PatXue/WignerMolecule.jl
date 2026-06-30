import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-rand"
tm.init_type = :rand
tm.init_T = 2.0
tm.bias = nothing
bias_type = Nothing

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ls = [24, 48]
Ts = 0.01:0.01:0.15
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 24)
    tm.thermalization = 75000 * div(L, 24)
    tm.binsize = 500 * div(L, 24)
    for T in Ts
        tm.T = T
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)