import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-high-t"
tm.init_type = :afm_fe
tm.bias = nothing

tm.sweeps = 20000
tm.thermalization = 20000
tm.binsize = 100

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
βs = 0.25:1.0:15.25
Ls = [8]
for (β, L) in Iterators.product(βs, Ls)
    tm.Lx = tm.Ly = L
    tm.β = β
    tm.T = 1/β
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)