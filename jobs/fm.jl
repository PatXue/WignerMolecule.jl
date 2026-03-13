import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "fm"

tm.sweeps = 50000
tm.thermalization = 50000
tm.binsize = 500
tm.init_type = :const

fm_bias(_, _) = [0, 0, 1]
tm.bias = fm_bias
bias_type = typeof(fm_bias)
tm.init_B = 20.0
tm.B = 0
JSON.lower(f::bias_type) = f(1, 1)

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 5, 20, 9)]...)
Ts = 0.5:0.5:7
Ls = [20, 40, 80]
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * L/20
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