import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "jp9"
tm.init_type = :rand
tm.init_T = 3.0
tm.B = 0.01

Ls = [24, 48]
Ts = [0.25, 0.5, 0.75, 1.1, 1.25, 1.5]
Jzs = 0.5:0.1:1.5
for (Jz, T, L) in Iterators.product(Jzs, Ts, Ls)
    tm.sweeps = 20000
    tm.thermalization = 20000
    tm.binsize = 200
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
