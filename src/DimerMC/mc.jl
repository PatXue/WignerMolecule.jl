struct DimerMC{AlgType} <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams
    fug::Float64        # Fugacity of a dimer

    spins::PeriodicMatrix{SVector{2,Int}}  # Matrix holding position (x,y) of entangled partner
    spinscopy::PeriodicMatrix{SVector{2,Int}}
    monospins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}
    monomers::BitSet

    sks::PeriodicArray{ComplexF64, 3}
    ηks::PeriodicArray{ComplexF64, 3}       # Fourier transformed ηs
    corr_rad::Int                   # Radius around which to sum correlations
end

function DimerMC(; T, init_T, wigparams, fug, Lx, Ly, algtype=:Heatbath, corr_rad=0)
    init_ss = fill(zeros(SVector{2,Int}), (Lx, Ly))
    copy_ss = algtype == :Worm ? copy(init_ss) : zeros(0, 0)
    init_ssmono = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))

    return DimerMC{algtype}(
        T, init_T, wigparams, fug,
        init_ss, copy_ss, init_ssmono, init_ηs, BitSet(1:(Lx*Ly)),
        Array{ComplexF64}(undef, (Lx, Ly, 4)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        corr_rad
    )
end

function DimerMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]
    fug = try
        params[:fug]
    catch
        exp(-params[:mu] / T)
    end
    algtype = get(params, :algtype, :Heatbath)
    corr_rad = get(params, :corr_rad, 0)

    return DimerMC(; T, init_T, wigparams, fug, Lx, Ly, algtype, corr_rad)
end
