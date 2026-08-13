import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "er7"
tm.algtype = :Heatbath
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

inits = [:stripe, :afm_fe, :afm_afe]
Ls = [48]
Ts = [0.5, 1.0, 2.0]
ams = 4:11
for (am, T, L, init) in Iterators.product(ams, Ts, Ls, inits)
    tm.init_type = init
    tm.Lx = tm.Ly = L
    tm.am = am
    tm.wigparams = WignerParams("all_params.jld2", 7, am, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)