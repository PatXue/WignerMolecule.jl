import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-bias"
tm.init_type = :afm_afe

afm_bias(x, y) = [0, 0, (-1)^(x + div(y,2))]
tm.bias = afm_bias
bias_type = typeof(afm_bias)
tm.init_B = 1.0
JSON.lower(f::bias_type) = f(1, 1)

raw_params = load_object("all_params.jld2")[(45, 11, 20, 7)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Bs = 0.0:0.005:0.05
Ts = [0.04, 0.045, 0.05]
Ls = [20, 40, 80]
for (B, T, L) in Iterators.product(Bs, Ts, Ls)
    tm.sweeps = 50000 * div(L, 20)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
    tm.Lx = tm.Ly = L
    tm.T = T
    tm.B = B
    task(tm)
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)