struct WignerParams
    J_SS::ComplexF64        # Spin-spin coupling
    J_EzEz_SS::ComplexF64   # Spin-ηz coupling
    J_EzEz::ComplexF64      # ηz coupling
    J_EAM_SS::ComplexF64    # Spin-η weird coupling (J+)
    J_EMEP_SS::ComplexF64   # Spin-η± coupling
    J_EMEM_SS::ComplexF64   # Spin-η- coupling
    J_EMEP::ComplexF64      # η± coupling
    J_EMEM::ComplexF64      # η- coupling
end

function WignerParams(paramfile, e_r, a_M, ϕ=45, d_g=20; H0=nothing)
    WignerParams(paramfile, (ϕ, e_r, d_g, a_M); H0)
end

function WignerParams(paramfile, idx; H0=nothing)
    raw_params = load_object(paramfile)[idx]
    params = [raw_params[i] for i in 1:8]
    H0 = isnothing(H0) ? norm(params) : H0
    norm_params = params ./ H0
    return WignerParams(norm_params...)
end

# Default WignerParams values (for testing)
const default_params = WignerParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

norm2(v) = sum(abs2.(v))
function LinearAlgebra.norm(wp::WignerParams)
    sqrt(norm2([wp.J_SS, wp.J_EAM_SS, wp.J_EMEM_SS, wp.J_EMEP_SS, wp.J_EzEz_SS, wp.J_EMEM, wp.J_EMEP, wp.J_EzEz]))
end

const SpinVector = SVector{3, Float64}
# Note: Using temperature in units of energy (k_B = 1)
mutable struct WignerMC{AlgType, BiasType} <: AbstractMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams
    H0::Float64         # Energy scale of params
    B::Float64          # Bias field magnitude
    init_B::Float64     # Thermalization bias field
    bias::BiasType      # Bias field pattern (normalized)

    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}

    spinks::PeriodicArray{ComplexF64, 3}    # Fourier transformed spins
    ηks::PeriodicArray{ComplexF64, 3}       # Fourier transformed ηs
    corr_rad::Int                           # Radius around which to sum correlations

    f::Float64              # Each sweep flips f*N spins
    Nc::Ref{Expectation}    # Avg cluster size
end

function WignerMC(; T=1.0, init_T=1.0, wigparams=default_params, B=0.0, init_B=0.0,
    bias=nothing, Lx=48, Ly=48, corr_rad=0, algtype=:Metropolis, f=1)
    biastype = typeof(bias)
    init_spins = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return WignerMC{algtype, biastype}(
        T, init_T, wigparams, norm(wigparams), B, init_B, bias, init_spins, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        corr_rad, f, Ref{Expectation}(Expectation(0.0, 0))
    )
end

function WignerMC(params::AbstractDict)
    algtype = get(params, :algtype, :Metropolis)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]
    corr_rad = get(params, :corr_rad, 0)
    f = get(params, :f, 1)

    B = get(params, :B, 0.0)
    init_B = get(params, :init_B, B)
    bias = get(params, :bias, nothing)

    return WignerMC(; T, init_T, wigparams, B, init_B, bias, Lx, Ly, corr_rad, algtype, f)
end
