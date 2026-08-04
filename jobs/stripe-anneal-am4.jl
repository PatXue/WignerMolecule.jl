import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using WignerMolecule

tm = TaskMaker()
jobname = "stripe-anneal-am4"

tm.init_type = :stripe
tm.sweeps = 100000
tm.thermalization = 100000
tm.binsize = 500

stripe_bias(x, _) = [0, 0, (-1)^(div(x, 2))]
tm.bias = stripe_bias
bias_type = typeof(stripe_bias)
tm.B = 0.0
tm.init_B = 5.0
JSON.lower(f::bias_type) = f(1, 1)

Ls = [48]
ers = 5:11
Ts = 0.3:0.3:4.5
for (T, er, L) in Iterators.product(Ts, ers, Ls)
    tm.Lx = tm.Ly = L
    tm.er = er
    tm.wigparams = WignerParams("all_params.jld2", er, 6, H0=1)
    tm.T = T
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)