import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-param"

L = 40
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :stripe

param_d = load_object("all_params.jld2")
tm.T = 0.1
params = Iterators.flatten((
    [(4, e) for e in 5:11],
    [(5, e) for e in 5:8],
    [(6, e) for e in 5:8],
    [(7, e) for e in 5:7],
    [(8, 5), (9, 5)]
),)
for param in params
    tm.a_m = param[1]
    tm.epsilon = param[2]
    raw_params = collect(param_d[(45, param[2], 20, param[1])])
    for i in (1, 2, 3, 6, 8)
        raw_params[i] = real(raw_params[i])
    end
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