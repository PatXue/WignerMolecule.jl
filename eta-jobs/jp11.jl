import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp11"
tm.init_type = :stripe
tm.B = 0.01
tm.allchis = true

Ls = [24]
Jzs = [0.5, 0.95, 1.0, 1.05, 1.125, 1.2]
Ts = 0.3:0.1:1.7
for (T, Jz, L) in Iterators.product(Ts, Jzs, Ls)
    tm.sweeps = 40000
    tm.thermalization = 40000
    tm.binsize = 400
    tm.wigparams = EtaParams(Jz, 1.1)
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
