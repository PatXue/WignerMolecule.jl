const J = 1     # In-plane coupling, set to 1 by convention
struct EtaParams
    Jzz::Float64    # z-component coupling
    Jp::Float64     # Frustrating coupling
end

struct EtaMC <: AbstractMC
    T::Float64          # Temperature
    params::EtaParams

    Ss::PeriodicMatrix{SpinVector}  # Spins
    Sks::Array{ComplexF64, 3}       # Fourier transformed spins
end

