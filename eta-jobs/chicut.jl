import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "chicut"
tm.B = 0.001
tm.allchis = true

Ls = [60]
Ts = [0.1, 0.75, 1.25]
Jps = Iterators.flatten((range(0.0, 0.9, 5), range(1.1, 2.0, 5)))
diffs = [-0.15, -0.05, 0.05, 0.1, 0.2]
for (d, Jp, T, L) in Iterators.product(diffs, Jps, Ts, Ls)
    tm.sweeps = 50000
    tm.thermalization = 50000
    tm.binsize = 250

    tm.wigparams = EtaParams(Jz, 0.9)
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jp = Jp
    tm.Jz = max(1.0, Jp) + d
    tm.init_type = (d > 0) ? :fm : ((Jp > 1) ? :stripe : :fe)
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
