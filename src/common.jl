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
    elseif init_type == :afm_afe
        init_afm_afe!(mc.spins, mc.ηs)
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
    E = total_energy(mc, mc.B) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in corr_posns
        pos = f(Lx, Ly)
        s = mc.spinks[pos..., :]
        η = mc.ηks[pos..., :]
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("ηk_", f), η)
        measure!(ctx, Symbol("sk_corr_", f), norm2(s))
        measure!(ctx, Symbol("ηk_corr_", f), η*η')

        measure!(ctx, Symbol("ηk_corr_", f, "_a1"), abs2(η[1]))
        measure!(ctx, Symbol("ηk_corr_", f, "_a2"), abs2(-η[1]/2 + η[2]*√3/2))
        measure!(ctx, Symbol("ηk_corr_", f, "_a3"), abs2(-η[1]/2 - η[2]*√3/2))
    end

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC{AlgType, BiasType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType, BiasType}
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N / T * (mag2 - mag^2)
    end

    for f in corr_posns
        evaluate!(eval, Symbol("χs_", f), (Symbol("sk_", f), Symbol("sk_corr_", f))) do sk, sk2
            N / T * (sk2 - sum(abs2.(sk)))
        end
        evaluate!(eval, Symbol("χηx_", f), (Symbol("ηk_", f), Symbol("ηk_corr_", f))) do ηk, ηk2
            N / T * real(ηk2[1,1] - abs2(ηk[1]))
        end
        evaluate!(eval, Symbol("χηy_", f), (Symbol("ηk_", f), Symbol("ηk_corr_", f))) do ηk, ηk2
            N / T * real(ηk2[2,2] - abs2(ηk[2]))
        end
        evaluate!(eval, Symbol("χηz_", f), (Symbol("ηk_", f), Symbol("ηk_corr_", f))) do ηk, ηk2
            N / T * real(ηk2[3,3] - abs2(ηk[3]))
        end
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

function Carlo.parallel_tempering_log_weight_ratio(mc, param_name::Symbol, new_val)
    if param_name != :T
        throw(ArgumentError("Parallel tempering only supported for T"))
    end
    return -(1/new_val - 1/mc.T) * total_energy(mc)
end
function Carlo.parallel_tempering_change_parameter!(mc, param_name::Symbol, new_val)
    if param_name != :T
        throw(ArgumentError("Parallel tempering only supported for T"))
    end
    mc.T = new_val
end
