import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jz2"
tm.init_type = :rand
tm.B = 0.001

tm.sweeps = 200000
tm.thermalization = 200000
tm.binsize = 500

Ls = [72]
Ts = Base.LogRange(1e-4, 1e-1, 5)
Jps = range(1.5, 2, 20)
for (Jp, T, L) in Iterators.product(Jps, Ts, Ls)
    tm.wigparams = EtaParams(2, Jp)
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jp = Jp
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
