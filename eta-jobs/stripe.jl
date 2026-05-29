import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "stripe"
tm.wigparams = EtaParams(0.5, 2)
tm.init_type = :stripe

Ls = [20]
Ts = 0.1:0.1:1.0
for (T, L) in Iterators.product(Ts, Ls)
    tm.sweeps = 20000 * div(L, 20)
    tm.thermalization = 20000 * div(L, 20)
    tm.binsize = 100 * div(L, 20)
    tm.T = T
    tm.Lx = tm.Ly = L
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
