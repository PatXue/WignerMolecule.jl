import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp15"
tm.B = 0.01

Ls = [24]
Jzs = [1.0, 1.45, 1.55, 2.0]
Ts = 0.3:0.1:1.7
for (T, Jz, L) in Iterators.product(Ts, Jzs, Ls)
    tm.sweeps = 40000
    tm.thermalization = 40000
    tm.binsize = 400
    tm.init_type = (Jz > 1.5) ? :fm : :stripe

    tm.wigparams = EtaParams(Jz, 1.5)
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
