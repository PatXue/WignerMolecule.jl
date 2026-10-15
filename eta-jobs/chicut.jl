import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "chicut"
tm.B = 0.001
tm.allchis = true
tm.corr_rad = 2
tm.init_type = :rand
tm.init_T = 5.0
tm.annealT = 0.0

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 500

Ls = [60]
Ts = [0.1, 0.5, 0.75, 1.0, 1.25]
Jzs = [0.4, 0.8, 0.9, 0.95, 0.975, 1.025, 1.05, 1.1, 1.5, 2.0]
for (Jz, T, L) in Iterators.product(Jzs, Ts, Ls)
    if Jz <= 1.0
        Jps = Iterators.flatten((range(0.1, 0.9, 5), 0.95, 1.05, range(1.1, 1.5, 3)))
    else
        Jps = Iterators.flatten((range(0.1, Jz-0.1, 5), Jz-0.05, Jz+0.05, range(Jz+0.1, Jz+0.5, 3)))
    end
    tm.T = T
    tm.Lx = tm.Ly = L
    tm.Jz = Jz
    for Jp in Jps
        tm.Jp = Jp
        tm.wigparams = EtaParams(Jz, Jp)
        task(tm)
    end
end

job = JobInfo("$jobname", EtaMC;
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)
