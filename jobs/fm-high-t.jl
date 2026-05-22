import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "fm"

tm.sweeps = 20000
tm.thermalization = 20000
tm.binsize = 100
tm.init_type = :const

fm_bias(_, _) = [0, 0, 1]
tm.bias = fm_bias
bias_type = typeof(fm_bias)
tm.init_B = 20.0
tm.B = 0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 5, 20, 9)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.6:0.1:2.5
Ls = [8]
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