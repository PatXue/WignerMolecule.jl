import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "fm-bias"

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000
tm.init_type = :rand

fm_bias(B) = (_, _) -> [0, 0, B]
tm.bias = fm_bias(1.0)
bias_type = typeof(tm.bias)

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 5, 20, 9)]...)
tm.T = 0.5
tm.Lx = tm.Ly = 80
Bs = 0.0:1.0:10.0
for B in Bs
    tm.B = B
    tm.bias = fm_bias(B)
    spins_dir = "$jobname.data/$(current_task_name(tm))"
    tm.outdir = spins_dir
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)