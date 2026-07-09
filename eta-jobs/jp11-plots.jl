import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp11-plots"
tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.B = 0.001

Ls = [24]
Jzs = [0.9, 1.05, 1.2, 1.4]
Ts = 0.1:0.1:2.0
for (T, Jz, L) in Iterators.product(Ts, Jzs, Ls)
    tm.init_type = (Jz > 1.1) ? :fm : :stripe
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
