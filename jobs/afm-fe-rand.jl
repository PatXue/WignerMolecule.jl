import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-rand"
tm.init_type = :rand
tm.bias = nothing
bias_type = Nothing

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.005:0.005:0.1
Ls = [20, 40]
for (T, L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 20)
    tm.thermalization = 50000 * div(L, 20)
    tm.binsize = div(50000 * div(L,20), 100)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)