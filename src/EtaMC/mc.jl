const J = 1     # In-plane coupling, set to 1 by convention
struct EtaParams
    Jzz::Float64    # z-component coupling
    Jp::Float64     # Frustrating coupling
end

struct EtaMC <: AbstractMC
    T::Float64
    params::EtaParams
    B::Float64
    spins::PeriodicMatrix{SpinVector}
    spinks::Array{ComplexF64, 3}        # Fourier transformed spins
    chis::Matrix{ComplexF64}
end

function EtaMC(; T=1.0, B=0.0, wigparams=EtaParams(0.0,0.0), Lx=24, Ly=24)
    init = fill(zeros(SpinVector), (Lx, Ly))
    return EtaMC(
        T, wigparams, B, init,
        Array{ComplexF64, 3}(undef, (Lx, Ly, 3)),
        Matrix{ComplexF64}(undef, Lx, Ly)
    )
end

function EtaMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    B = get(params, :B, 0.0)
    return EtaMC(; T, B, Lx, Ly, wigparams)
end
