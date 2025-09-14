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
const default_params = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

# Note: Using temperature in units of energy (k_B = 1)
struct WignerMC{AlgType} <: AbstractMC
    T::Float64     # Temperature
    params::WignerParams
    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}
end

function WignerMC{AlgType}(; T=1.0, wigparams=default_params, Lx=40, Ly=40) where {AlgType}
    return WignerMC{AlgType}(T, WignerParams(wigparams...),
                             fill(zeros(SpinVector), (Lx, Ly)))
end

function WignerMC{AlgType}(params::AbstractDict) where {AlgType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    return WignerMC{AlgType}(; T, wigparams, Lx, Ly)
end

function Carlo.init!(mc::WignerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = params[:init_type]
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
        end
    else
        rand!(ctx.rng, mc.spins)
        rand!(ctx.rng, mc.ηs)
    end
    return nothing
end

function Carlo.sweep!(mc::WignerMC, ctx::Carlo.MCContext)
    Carlo.sweep!(mc, ctx.rng)
end

const ω::ComplexF64 = exp(im * 2π/3)
# Calculate the energy at a lattice site (x, y) if it had spin s and
# pseudospin η
function energy(mc::WignerMC, s::SpinVector, η::SpinVector, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1), (x-1, y), (x-1, y+1), (x, y+1))
    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]]
        ηj = mc.ηs[nns[j]]

        # η raising and lowering operators
        η_m = η[1] + 1.0im*η[2]
        ηj_p = ηj[1] - 1.0im*ηj[2]
        ηj_m = ηj[1] + 1.0im*ηj[2]

        E_spin = 0.0 + 0.0im
        E_η = 0.0 + 0.0im

        # η-only energy
        E_η +=   J_EzEz *     η[3] * ηj[3]
        E_η += 2*J_EMEP *     η_m * ηj_p
        E_η += 2*J_EMEM * ν * η_m * ηj_m

        # η-S energy
        E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
        E_spin += 2*J_EMEP_SS *     η_m * ηj_p
        E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
        E_spin += J_SS
        E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)

        E += real(E_spin * (s⋅sj) + E_η)
    end

    return E
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::WignerMC, x, y)
    s = mc.spins[x, y]
    nn = mc.spins[x+1, y] + mc.spins[x, y+1]
    nnna = mc.spins[x+1, y+1]
    nnnb = mc.spins[x+1, y-1]
    H0 = s ⋅ (mc.J1 * nn + mc.J2a * nnna + mc.J2b * nnnb)
    biquad = mc.K * (s ⋅ mc.spins[x+1, y])^2
    return H0 + biquad
end

function Carlo.measure!(mc::WignerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Magnetization per lattice site
    mag = norm(sum(mc.spins)) / N
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Mag4, mag^4)

    # Energy per lattice site
    energy = 0.0
    for y in 1:Ly
        for x in 1:Lx
            energy += half_energy(mc, x, y)
        end
    end
    energy /= N
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC{AlgType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType}
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N / T * (mag2 - mag^2)
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end

function Carlo.write_checkpoint(mc::WignerMC, out::HDF5.Group)
    out["spins"] = mc.spins
    out["etas"] = mc.ηs
    return nothing
end
function Carlo.read_checkpoint!(mc::WignerMC, in::HDF5.Group)
    raw_spins = read(in, "spins")
    raw_ηs = read(in, "etas")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_spins)
    mc.ηs .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_ηs)
    return nothing
end
