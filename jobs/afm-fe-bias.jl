import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using LinearAlgebra
using JLD2
using WignerMolecule

tm = TaskMaker()
jobname = "afm-fe-bias"

tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 1000

afm_bias(B) = (x, _) -> [0, 0, B * (-1)^x]
tm.bias = afm_bias(1.0)
bias_type = typeof(tm.bias)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 10)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
tm.init_T = 10
tm.T = 0.01
Ls = [20]
Bs = 0.0:2.5:40.0
for L in Ls
    tm.Lx = tm.Ly = L
    for B in Bs
        spins_dir = "$jobname.data/$(current_task_name(tm))"
        tm.outdir = spins_dir
        tm.B = B
        tm.bias = afm_bias(B)
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)