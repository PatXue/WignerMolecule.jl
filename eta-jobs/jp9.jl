import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp9"
tm.init_type = :fe
tm.B = 0.001

Ls = [24, 48]
Jzs = [0.75, 0.9, 0.95, 1.05, 1.2, 1.5]
Ts = 0.3:0.1:1.7
for (T, Jz, L) in Iterators.product(Ts, Jzs, Ls)
    tm.sweeps = 40000
    tm.thermalization = 40000
    tm.binsize = 400
    tm.init_type = (Jz > 1) ? :fm : :fe

    tm.wigparams = EtaParams(Jz, 0.9)
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
