import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-crit"
tm.init_type = :stripe
tm.algtype = :Cluster
tm.sweeps = 200000
tm.thermalization = 400000
tm.binsize = 2000

n = 20
tm.wigparams = WignerParams("all_params.jld2", 5, 6)
Ls = [12, 16, 20, 24, 36]
Tbounds = [(0.081, 0.084), (0.0805, 0.082), (0.0805, 0.082), (0.08, 0.0815), (0.0802, 0.081)]
for (L, (Tmin, Tmax)) in Iterators.zip(Ls, Tbounds)
    Ts = range(Tmin, Tmax, n)
    tm.parallel_tempering = (
        mc = WignerMC,
        parameter = :T,
        values = Ts,
        interval = 1
    )
    tm.Lx = tm.Ly = L
    task(tm)
end

job = JobInfo("$jobname", ParallelTemperingMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = n
)
start(job, ARGS)