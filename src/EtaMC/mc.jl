const J = 1     # In-plane coupling, set to 1 by convention
struct EtaParams
    Jzz::Float64    # z-component coupling
    Jp::Float64     # Frustrating coupling
end

struct EtaMC <: AbstractMC
    T::Float64          # Temperature
    wigparams::EtaParams

    Ss::PeriodicMatrix{SpinVector}  # Spins
    Sks::Array{ComplexF64, 3}       # Fourier transformed spins
end

function EtaMC(; T=1.0, wigparams, Lx=20, Ly=20)
    init = fill(zeros(SpinVector), (Lx, Ly))
    return EtaMC(T, wigparams, init, Array{ComplexF64, 3}(undef, (Lx, Ly, 3)))
end

function EtaMC(params::AbstractDict)
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    return EtaMC(; T, Lx, Ly, wigparams)
end
