import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "fecut"
tm.B = 0.001
tm.allchis = true
tm.corr_rad = 2
tm.init_type = :rand
tm.init_T = 5.0
tm.annealT = 0.0

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

Ls = [24]
Ts = [0.1, 0.5, 0.75, 1.0, 1.25]
diffs = [-0.15, -0.1, -0.05, -0.025, 0.025, 0.05, 0.1, 0.15]
Jps = range(0.1, 0.95, 10)
for (Jp, d, T, L) in Iterators.product(Jps, diffs, Ts, Ls)
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jp = Jp
    tm.diff = d
    tm.Jz = max(1.0, Jp) + d
    tm.wigparams = EtaParams(tm.Jz, Jp)
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
