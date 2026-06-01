import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "gap"
tm.wigparams = EtaParams(0.5, 0.5)

Ls = [20]
Ts = [0.1]
Jzs = 0.0:0.1:2
for (Jz, T, L) in Iterators.product(Jzs, Ts, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 100
    tm.wigparams = EtaParams(Jz, 1.0)
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jz = Jz
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
