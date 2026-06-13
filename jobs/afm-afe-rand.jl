import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-rand"
tm.init_type = :rand
tm.init_T = 0.5
tm.bias = nothing
bias_type = Nothing

raw_params = load_object("all_params.jld2")[(45, 11, 20, 7)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.0025:0.0025:0.05
Ls = [24, 48]
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000 * div(L, 24)
    tm.thermalization = tm.sweeps
    tm.binsize = div(tm.sweeps, 100)
    for T in Ts
        tm.T = T
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, bias_type};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)