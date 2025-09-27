import Pkg
Pkg.activate("..")

using WignerMolecule
using Carlo
using Carlo.JobTools

tm = TaskMaker()

L = 20
tm.Lx = tm.Ly = L
tm.sweeps = 20000
tm.thermalization = 0
tm.binsize = 100
tm.init_type = :afm_fe

tm.savefreq = 5000

tm.wigparams = WignerParams(
    12.006059575026349,
    -2.655552440850471,
    -0.626818858667411,
    -7.583065735252017+7.793472257363657im,
    0.4778039991600792-5.6720434985650225im,
    6.241388905807156,
    0.3787671354263329-1.3128412213891711im,
    1.8423139895402252,
    0
)
Ts = 0.05:0.05:1.0
for T in Ts
    tm.T = T
    spins_dir = "afm-fe.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("afm-fe", WignerMC{:Metropolis};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)