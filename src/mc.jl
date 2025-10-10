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
    H::Float64              # Bias field
end

# Default WignerParams values (for testing)
const default_params = WignerParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-3)

const SpinVector = SVector{3, Float64}
# Note: Using temperature in units of energy (k_B = 1)
struct WignerMC{AlgType} <: AbstractMC
    T::Float64     # Temperature
    params::WignerParams
    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}

    spinks::Array{ComplexF64, 3}    # Fourier transformed spins
    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs
    spink_corrs::Matrix{Float64}    # Correlations of Fourier spins
    ηk_corrs::Array{ComplexF64, 4}  # Correlations of Fourier ηs

    outdir::String # Output directory for local spin current plots
    savefreq::Int  # No. of sweeps between saving local spin current
end

function WignerMC{AlgType}(; T=1.0, wigparams=default_params, Lx=40, Ly=40,
    outdir=".", savefreq=0) where {AlgType}
    init_spins = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return WignerMC{AlgType}(
        T, wigparams, init_spins, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Matrix{Float64}(undef, (Lx, Ly)),
        Array{ComplexF64}(undef, (Lx, Ly, 3, 3)),
        outdir, savefreq
    )
end

function WignerMC{AlgType}(params::AbstractDict) where {AlgType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    outdir = get(params, :outdir, ".")
    savefreq = get(params, :savefreq, 0)
    return WignerMC{AlgType}(; T, wigparams, Lx, Ly, outdir, savefreq)
end
