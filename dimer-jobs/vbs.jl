import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "vbs"

raw_params = load_object("all_params.jld2")[(45, 11, 20, 6)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)

Ts = 0.01:0.01:0.1
Ls = [24]
for (T, L) in Iterators.product(Ts, Ls)
    tm.sweeps = 25000 * div(L, 24)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
    tm.Lx = tm.Ly = L
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)