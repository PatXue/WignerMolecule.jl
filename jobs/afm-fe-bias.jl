import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-bias"

tm.sweeps = 20000
tm.thermalization = 20000
tm.binsize = 200

afm_bias(B) = (x, _) -> [0, 0, B * (-1)^x]
tm.bias = afm_bias(1.0)
bias_type = typeof(tm.bias)

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 11, 20, 10)]...)
tm.init_T = 10
Ts = 0.0:0.1:1.0
Ls = [20]
for L in Ls
    tm.Lx = tm.Ly = L
    for T in Ts
        tm.T = max(T, 0.01)
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)