import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "am9"
tm.algtype = :Heatbath
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500

inits = [:fm, :afm_fe]
Ls = [48]
ers = 5:11
Ts = 0.2:0.2:2.0
for (T, er, L, init) in Iterators.product(Ts, ers, Ls, inits)
    tm.init_type = init
    tm.Lx = tm.Ly = L
    tm.er = er
    tm.wigparams = WignerParams("all_params.jld2", er, 9, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)