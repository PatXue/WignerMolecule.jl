import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jz_cut"
tm.init_type = :rand

Ls = [20]
Ts = [0.1]
Jps = 0.5:0.1:1.5
for (Jp, T, L) in Iterators.product(Jps, Ts, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 100
    tm.wigparams = EtaParams(0.4, Jp)
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
