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

function WignerParams(paramfile, e_r, a_M, ϕ=45, d_g=20; H0=nothing)
    raw_params = load_object(paramfile)[(ϕ, e_r, d_g, a_M)]
    params = [raw_params[i] for i in 1:8]
    for i in (1, 2, 3, 6, 8)
        try
            convert(Float64, params[i])
        catch e
            if isapprox(imag(params[i]), 0, atol=1e-5)
                params[i] = real(params[i])
            else
                throw(e)
            end
        end
    end
    H0 = isnothing(H0) ? norm(params) : H0
    norm_params = params ./ H0
    return WignerParams(norm_params...)
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

    spinks::PeriodicArray{ComplexF64, 3}    # Fourier transformed spins
    ηks::PeriodicArray{ComplexF64, 3}       # Fourier transformed ηs
    corr_rad::Int                           # Radius around which to sum correlations
end

function WignerMC(; T=1.0, init_T=1.0, wigparams=default_params, B=0.0, init_B=0.0,
    bias=nothing, Lx=48, Ly=48, corr_rad=0, algtype=:Metropolis)
    biastype = typeof(bias)
    init_spins = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return WignerMC{algtype, biastype}(
        T, init_T, wigparams, B, init_B, bias, init_spins, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        corr_rad
    )
end

function WignerMC(params::AbstractDict)
    algtype = get(params, :algtype, :Metropolis)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    init_T = get(params, :init_T, T)
    wigparams = params[:wigparams]
    corr_rad = get(params, :corr_rad, 0)

    B = get(params, :B, 0.0)
    init_B = get(params, :init_B, B)
    bias = get(params, :bias, nothing)

    return WignerMC(; T, init_T, wigparams, B, init_B, bias, Lx, Ly, corr_rad, algtype)
end
