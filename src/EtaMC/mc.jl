const J = 1     # In-plane coupling, set to 1 by convention
struct EtaParams
    Jzz::Float64    # z-component coupling
    Jp::Float64     # Frustrating coupling
end

struct EtaMC <: AbstractMC
    T::Float64
    params::EtaParams
    B::Float64
    init_T::Float64

    spins::PeriodicMatrix{SpinVector}
    spinks::Array{ComplexF64, 3}        # Fourier transformed spins

    chis::Matrix{ComplexF64}
    allchis::Bool
end

function EtaMC(; T=1.0, init_T=1.0, B=0.0, wigparams=EtaParams(0.0,0.0), Lx=24, Ly=24, allchis=false)
    init = fill(zeros(SpinVector), (Lx, Ly))
    return EtaMC(
        T, wigparams, B, init_T, init,
        Array{ComplexF64, 3}(undef, (Lx, Ly, 3)),
        Matrix{ComplexF64}(undef, Lx, Ly),
        allchis
    )
end

function EtaMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    B = get(params, :B, 0.0)
    init_T = get(params, :init_T, T)
    allchis = get(params, :allchis, false)
    return EtaMC(; T, init_T, B, Lx, Ly, wigparams, allchis)
end
