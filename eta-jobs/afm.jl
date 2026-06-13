import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "afm"
tm.init_type = :afm
tm.wigparams = EtaParams(-10,0)

Ls = [24]
Ts = 0.2:0.2:2.0
for (T, L) in Iterators.product(Ts, Ls)
    tm.sweeps = 50000
    tm.thermalization = 50000
    tm.binsize = 250
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
