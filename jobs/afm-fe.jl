import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools

tm = TaskMaker()

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 40000
tm.thermalization = 60000
tm.binsize = 100
tm.init_type = :afm_fe

tm.wigparams = WignerParams(
    21.7432,
    -7.96475,
    -2.15716,
    (-13.5052+11.9152im),
    (1.95998-8.79458im),
    11.1589,
    (0.778735-1.96819im),
    3.55356,
)
Ts = 0.0:0.025:0.3
for T in Ts
    # tm.thermalization = T < 0.3 ? 40000 : 20000
    tm.T = max(T, 0.01)
    spins_dir = "afm-fe-center.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("afm-fe-center", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)