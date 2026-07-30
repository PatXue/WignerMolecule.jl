import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "vbs-inf"
tm.init_type = :rand
tm.wigparams = WignerParams("all_params.jld2", 10, 6)

tm.sweeps = 25000
tm.thermalization = 25000
tm.binsize = 250
tm.T = Inf
fugs = 0.0:0.1:1.9
Ls = [24, 48]
for (fug, L) in Iterators.product(fugs, Ls)
    tm.Lx = tm.Ly = L
    tm.fug = fug
    task(tm)
end

job = JobInfo("$jobname", DimerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)