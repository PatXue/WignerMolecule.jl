struct EtaParams
    J::ComplexF64      # in-plane coupling
    Jzz::Float64    # z-component coupling
    Jp::Float64     # Frustrating coupling
end

EtaParams(Jzz, Jp) = EtaParams(1, Jzz, Jp)
function EtaParams(params::WignerParams)
    J_EzEz_SS = params.J_EzEz_SS
    J_EMEP_SS = params.J_EMEP_SS
    J_EMEM_SS = params.J_EMEM_SS
    J_EzEz = params.J_EzEz
    J_EMEP = params.J_EMEP
    J_EMEM = params.J_EMEM
    J = (J_EMEP_SS/4 + J_EMEP) / 4
    Jzz = (J_EzEz_SS/4 + J_EzEz) / 4
    Jp = (J_EMEM_SS/4 + J_EMEM) / 4
    return EtaParams(J, Jzz, Jp)
end

struct EtaMC{AlgType} <: AbstractMC
    T::Float64
    params::EtaParams
    B::Float64
    init_T::Float64

    spins::PeriodicMatrix{SpinVector}
    spinks::PeriodicArray{ComplexF64, 3}        # Fourier transformed spins

    chis::Matrix{ComplexF64}
    corr_rad::Int
end

function EtaMC(; T=1.0, init_T=1.0, B=0.0, etaparams=EtaParams(0.0,0.0), algtype=:Metropolis,
    Lx=24, Ly=24, corr_rad=0)
    init = fill(zeros(SpinVector), (Lx, Ly))
    return EtaMC{algtype}(
        T, etaparams, B, init_T, init,
        Array{ComplexF64, 3}(undef, (Lx, Ly, 3)),
        Matrix{ComplexF64}(undef, Lx, Ly),
        corr_rad
    )
end

function EtaMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    etaparams = params[:wigparams]

    algtype = get(params, :algtype, :Metropolis)
    B = get(params, :B, 0.0)
    init_T = get(params, :init_T, T)
    corr_rad = get(params, :corr_rad, 0)

    return EtaMC(; T, init_T, B, algtype, Lx, Ly, etaparams, corr_rad)
end
