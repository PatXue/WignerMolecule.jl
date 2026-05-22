# Direction of an entanglement at a site (a1 = xhat, proceeds counterclockwise)
@enum Bond a1 a2 a3 a4 a5 a6

struct DimerMC <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams

    spins::PeriodicMatrix{Bond}
    ηs::PeriodicMatrix{SpinVector}

    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs

    outdir::String # Output directory for spin plots
    savefreq::Int  # No. of sweeps between saving spins
end

function DimerMC(; T, init_T, wigparams, Lx, Ly, outdir="", savefreq=0)
    init_ss = fill(a1, (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return DimerMC(
        T, init_T, wigparams, init_ss, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
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
