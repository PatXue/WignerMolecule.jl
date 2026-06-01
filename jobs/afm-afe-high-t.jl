import Pkg
Pkg.activate("..")

using Carlo
using Carlo.JobTools
using JLD2
using JSON
using LinearAlgebra
using WignerMolecule

tm = TaskMaker()
jobname = "afm-afe-high-t"
tm.init_type = :afm_afe
tm.bias = nothing

raw_params = load_object("all_params.jld2")[(45, 11, 20, 7)]
norm_params = raw_params ./ norm(raw_params)
tm.wigparams = WignerParams(norm_params...)
Ts = 0.3:0.1:2.2
Ls = [12]
for L in Ls
    tm.Lx = tm.Ly = L
    tm.sweeps = 50000
    tm.thermalization = tm.sweeps
    tm.binsize = 100
    for T in Ts
        tm.T = max(T, 0.01)
        task(tm)
    end
end

job = JobInfo("$jobname", WignerMC{:Metropolis, Nothing};
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
)
start(job, ARGS)