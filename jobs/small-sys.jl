import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools

tm = TaskMaker()

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 50000
tm.thermalization = 0
tm.binsize = 100

tm.savefreq = 5000

Ts = (0.25, 4.0)
for T in Iterators.product(Ts)
    tm.T = T
    spins_dir = "small-sys.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("small-sys", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)