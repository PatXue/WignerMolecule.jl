import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp11"
tm.init_type = :stripe
tm.B = 0.01

Ls = [24, 48]
Jzs = [0.5, 1.0, 1.075, 1.125, 1.2, 1.5]
Ts = 0.3:0.1:1.5
for (T, Jz, L) in Iterators.product(Ts, Jzs, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 200
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
