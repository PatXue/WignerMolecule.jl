import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-coupl"

L = 40
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :stripe

param_d = load_object("all_params.jld2")
tm.T = 0.1
mults = 0.0:0.1:2.0
for mult in mults
    tm.mult = mult
    raw_params = collect(param_d[(45, 5, 20, 7)])
    raw_params[4] *= mult
    norm_params = raw_params ./ norm(raw_params)
    tm.wigparams = WignerParams(norm_params...)

    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis_η, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)