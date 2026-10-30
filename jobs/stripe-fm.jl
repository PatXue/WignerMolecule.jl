import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-fm"
tm.algtype = :Heatbath
tm.init_type = :rand
tm.init_T = 10.0
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

allparams = load_object("params_interpolated_quarter.jld2")
Ls = [60]
Ts = [0.5, 1.0, 2.0]
ams = range(8, 9, 10)
for (am, T, L) in Iterators.product(ams, Ts, Ls)
    tm.corr_rad = 2
    tm.Lx = tm.Ly = L
    tm.am = am
    if floor(am) == am
        raw = allparams[(am, 5)]
    else
        raw = (ceil(am)-am) * allparams[(floor(am), 5)] + (am-floor(am)) * allparams[(ceil(am), 5)]
    end
    tm.wigparams = WignerParams(raw...)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)