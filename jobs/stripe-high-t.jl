import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-high-t"
tm.init_type = :stripe
tm.bias = nothing

raw_params = load_object("all_params.jld2")[(45, 5, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ls = [8]
Ts = Iterators.flatten((0.04:0.01:0.09, 0.1:0.1:1.1))
for (T,L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 100
    tm.T = T
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)