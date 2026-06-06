import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "fm-high-t"

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :const
tm.bias = nothing

raw_params = load_object("all_params.jld2")[(45, 5, 20, 9)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = Iterators.flatten((0.05:0.025:0.1, 0.2:0.1:2.0))
Ls = [8]
for L in Ls
    tm.Lx = tm.Ly = L
    for T in Ts
        tm.T = T
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)