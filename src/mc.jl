# All energy units are in units of meV
struct WignerParams
    J_SS::Float64           # Spin-spin coupling
    J_EzEz_SS::Float64      # Spin-ηz coupling
    J_EzEz::Float64         # ηz coupling
    J_EAM_SS::ComplexF64    # Spin-η weird coupling (J+)
    J_EMEP_SS::ComplexF64   # Spin-η± coupling
    J_EMEM_SS::Float64      # Spin-η- coupling
    J_EMEP::ComplexF64      # η± coupling
    J_EMEM::Float64         # η- coupling
end

# Default WignerParams values (for testing)
const default_params = WignerParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

const SpinVector = SVector{3, Float64}
# Note: Using temperature in units of energy (k_B = 1)
struct WignerMC{AlgType, BiasType} <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams
    B::Float64          # Bias field magnitude
    init_B::Float64     # Thermalization bias field
    bias::BiasType      # Bias field pattern (normalized)

    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}

    spinks::Array{ComplexF64, 3}    # Fourier transformed spins
    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs

    outdir::String # Output directory for local spin current plots
    savefreq::Int  # No. of sweeps between saving local spin current
end

function WignerMC{AlgType, BiasType}(; T=1.0, init_T=1.0, wigparams=default_params,
    B, init_B, bias, Lx=40, Ly=40, outdir=".", savefreq=0) where {AlgType, BiasType}
    init_spins = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return WignerMC{AlgType, BiasType}(
        T, init_T, wigparams, B, init_B, bias, init_spins, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        outdir, savefreq
    )
end

function WignerMC{AlgType, BiasType}(params::AbstractDict) where {AlgType, BiasType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]

    B = get(params, :B, 0.0)
    init_B = get(params, :init_B, B)
    bias = get(params, :bias, nothing)

    outdir = get(params, :outdir, ".")
    savefreq = get(params, :savefreq, 0)

    return WignerMC{AlgType, BiasType}(;
        T, init_T, wigparams, B, init_B, bias, Lx, Ly, outdir, savefreq)
end
