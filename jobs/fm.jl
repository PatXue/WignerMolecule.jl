import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools

tm = TaskMaker()

L = 40
tm.Lx = tm.Ly = L
tm.sweeps = 20000
tm.thermalization = 0
tm.binsize = 100
tm.init_type = :const

tm.savefreq = 5000

tm.wigparams = WignerParams(-2.12742, -7.37151, -2.60026, (-1.5492-3.67457im), (2.22406+3.06319im), 0.672249, (0.215811+0.332566im), 0.815589, 1e-2)
Ts = 1:0.25:5
for T in Ts
    tm.T = T
    spins_dir = "fm.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("fm", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)