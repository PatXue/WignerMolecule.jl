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
    # Averaged adjacent dot products
    Dx0 = Dy0 = 0.0
    Dxπ = Dyπ = 0.0
    # Spin current
    spin_curr = zeros(3)
    x_hat = SVector(1, 0, 0)

    for y in 1:Ly
        for x in 1:Lx
            energy += half_energy(mc, x, y)

            s = mc.spins[x, y]
            sx = mc.spins[x+1, y]
            sy = mc.spins[x, y+1]
            x_dot = s ⋅ sx
            y_dot = s ⋅ sy
            Dx0 += x_dot
            Dy0 += y_dot
            Dxπ += x_dot * (-1)^(x+y)
            Dyπ += y_dot * (-1)^(x+y)

            spin_curr += s × sx
        end
    end
    energy /= N
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)
    Dx0 /= N
    Dy0 /= N
    measure!(ctx, :Dx0, Dx0)
    measure!(ctx, :Dy0, Dy0)
    Dxπ /= N
    Dyπ /= N
    measure!(ctx, :Dxπ, abs(Dxπ))
    measure!(ctx, :Dyπ, abs(Dyπ))
    spin_curr /= N
    measure!(ctx, :J_s, norm(spin_curr))
    measure!(ctx, :P, norm(x_hat × spin_curr))

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC{AlgType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType}
    T = params[:T]
    J2a = params[:J2a]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N * J2a/T * (mag2 - mag^2)
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
