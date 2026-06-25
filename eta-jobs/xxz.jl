import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "xxz"
tm.init_type = :rand
tm.init_T = 4.0

Ls = [48]
Ts = [0.25, 0.75, 1.1, 1.25, 1.5]
Jzs = 0.75:0.05:1.25
for (Jz, T, L) in Iterators.product(Jzs, Ts, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 100
    tm.wigparams = EtaParams(Jz, 0.0)
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
