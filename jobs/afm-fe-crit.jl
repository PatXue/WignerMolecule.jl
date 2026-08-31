import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-crit"
tm.init_type = :afm_fe
tm.sweeps = 200000
tm.thermalization = 300000
tm.binsize = 2000
tm.f = 2

tm.wigparams = WignerParams("all_params.jld2", 11, 10)
Ts = Iterators.flatten((0.051:0.00005:0.0515, 0.052))
Ls = [48, 72]
for (T, L) in Iterators.product(Ts, Ls)
    tm.Lx = tm.Ly = L
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", ClusterWigMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)