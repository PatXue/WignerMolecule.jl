import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "full-diagram-broad"
tm.B = 0.001
tm.corr_rad = 2
tm.init_type = :rand
tm.init_T = 5.0
tm.annealT = 0.0

tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 500

Ls = [60]
Ts = [0.8]
Jzs = range(0, 2, 40)
Jps = range(0, 2, 40)
for (Jp, Jz, T, L) in Iterators.product(Jps, Jzs, Ts, Ls)
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jz = Jz
    tm.Jp = Jp
    tm.wigparams = EtaParams(Jz, Jp)
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
