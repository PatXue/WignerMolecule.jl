import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "fm-bias"
tm.init_type = :const

fm_bias(_, _) = [0, 0, 1]
tm.bias = fm_bias
bias_type = typeof(fm_bias)
tm.init_B = 10.0
JSON.lower(f::bias_type) = f(1, 1)

tm.wigparams = WignerParams(load_object("all_params.jld2")[(45, 5, 20, 9)]...)
tm.Lx = tm.Ly = 20
Ls = [20, 40, 80]
Ts = [4.5, 5.5, 6.5]
Bs = 0.0:1.0:10.0
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 20)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
    for T in Ts
        tm.T = T
        for B in Bs
            tm.B = B
            spins_dir = "$jobname.data/$(current_task_name(tm))"
            tm.outdir = spins_dir
            task(tm)
        end
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)