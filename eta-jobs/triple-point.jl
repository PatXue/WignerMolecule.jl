import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "triple-point"
tm.init_type = :fm
tm.wigparams = EtaParams(1.0, 1.0)

Ls = [24, 48, 96]
Ts = 0.1:0.1:2.0
for (T, L) in Iterators.product(Ts, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 100
    tm.T = T
    tm.Lx = tm.Ly = L
    task(tm)
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
