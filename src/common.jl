# Initialize the WignerMC spins
function Carlo.init!(mc::WignerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = get(params, :init_type, :rand)

    rand!(ctx.rng, mc.spins)
    rand!(ctx.rng, mc.ηs)
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
            mc.ηs[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :afm_fe
        init_afm_fe!(mc.spins, mc.ηs)
    elseif init_type == :afm_fe_s
        init_afm_fe_s!(mc.spins)
    elseif init_type == :afm_fe_eta
        init_afm_fe_eta!(mc.ηs)
    elseif init_type == :stripe
        init_stripe!(mc.spins, mc.ηs)
    elseif init_type == :stripe_s
        init_stripe_s!(mc.spins)
    elseif init_type == :stripe_eta
        init_stripe_eta!(mc.ηs)
    end

    update_fourier!(mc)
    return nothing
end

function Carlo.measure!(mc::WignerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Magnetization per lattice site
    mag_v = sum(mc.spins) ./ N
    mag = norm(mag_v)
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Magx, mag_v[1])
    measure!(ctx, :Magy, mag_v[2])
    measure!(ctx, :Magz, mag_v[3])

    η = sum(mc.ηs) ./ N
    measure!(ctx, :ηx, η[1])
    measure!(ctx, :ηy, η[2])
    measure!(ctx, :ηz, abs(η[3]))
    measure!(ctx, :ηxy, sqrt(η[1]^2 + η[2]^2))

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

    update_fourier!(mc)
    measure!(ctx, :spink_corrs, mc.spink_corrs)
    measure!(ctx, :etak_corrs, mc.ηk_corrs)

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC{AlgType, B}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType, B}
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
