# Direction of an entanglement at a site (a1 = xhat, proceeds counterclockwise)
@enum Bond a1=0 a2=1 a3=2 a4=3 a5=4 a6=5 none

struct DimerMC <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams

    spins::PeriodicMatrix{Bond}
    visited::PeriodicMatrix{Bool}
    ηs::PeriodicMatrix{SpinVector}

    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs

    outdir::String # Output directory for spin plots
    savefreq::Int  # No. of sweeps between saving spins
end

function DimerMC(; T, init_T, wigparams, Lx, Ly, outdir="", savefreq=0)
    init_ss = fill(none, (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return DimerMC(
        T, init_T, wigparams,
        init_ss, Matrix(false, Lx, Ly),
        init_ηs, Array{ComplexF64}(undef, (Lx, Ly, 3)),
        outdir, savefreq
    )
end

function DimerMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]

    outdir = get(params, :outdir, ".")
    savefreq = get(params, :savefreq, 0)

    return DimerMC(; T, init_T, wigparams, Lx, Ly, outdir, savefreq)
end
