mutable struct DimerMC <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams
    fug::Float64        # Fugacity of a dimer

    spins::PeriodicMatrix{SVector{2,Int}}  # Matrix holding position (x,y) of entangled partner
    monospins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}
    monomers::BitSet

    sks::Array{ComplexF64, 3}
    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs

    etaonly::Bool

    outdir::String # Output directory for spin plots
    savefreq::Int  # No. of sweeps between saving spins
end

function DimerMC(; T, init_T, wigparams, fug, Lx, Ly, etaonly=false, outdir="", savefreq=0)
    init_ss = fill(zeros(SVector{2,Int}), (Lx, Ly))
    init_ssmono = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return DimerMC(
        T, init_T, wigparams, fug,
        init_ss, init_ssmono, init_ηs, BitSet(1:(Lx*Ly)),
        Array{ComplexF64}(undef, (Lx, Ly, 4)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        etaonly, outdir, savefreq
    )
end

function DimerMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]
    fug = get(params, :fug, 1.0)

    etaonly = get(params, :etaonly, false)
    outdir = get(params, :outdir, ".")
    savefreq = get(params, :savefreq, 0)

    return DimerMC(; T, init_T, wigparams, fug, Lx, Ly, etaonly, outdir, savefreq)
end
