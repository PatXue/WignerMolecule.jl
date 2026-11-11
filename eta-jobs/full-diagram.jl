import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using WignerMolecule

tm = TaskMaker()
jobname = "full-diagram"
tm.B = 0.001
tm.allchis = true
tm.corr_rad = 2
tm.init_type = :rand
tm.init_T = 5.0
tm.annealT = 0.0

tm.sweeps = 100000
tm.thermalization = 200000
tm.binsize = 500

Ls = [60]
Ts = [0.1, 0.5, 0.75, 1.0, 1.5]
Jzs = range(0.8, 2, 20)
for (Jz, T, L) in Iterators.product(Jzs, Ts, Ls)
    if Jz == 1.0
        Jz += 0.01
    end
    b = max(1.0, Jz)
    Jps = Iterators.flatten((
        range(b-0.8, b-0.3, 3),
        b-0.25, b-0.2, b-0.15, b-0.1, b-0.05, b-0.025, b+0.025, b+0.05, b+0.1, b+0.15, b+0.2
    ))
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
