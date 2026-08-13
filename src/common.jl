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
    measure!(ctx, :Mag4, mag^4)

    η = sum(mc.ηs) ./ N
    measure!(ctx, :ηz, abs(η[3]))
    measure!(ctx, :ηxy, sqrt(η[1]^2 + η[2]^2))

    # Energy per lattice site
    E = total_energy(mc) / N
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
        if mc.corr_rad != 0
            r = mc.corr_rad
            x, y = pos[1], pos[2]
            scorr = sum(abs2, mc.spinks[x-r:x+r, y-r:y+r, :])
            ηcorr = sum(ηk -> ηk*ηk', eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
            measure!(ctx, Symbol("sk_corr_near_", f), scorr)
            measure!(ctx, Symbol("ηk_corr_near_", f), ηcorr)
        end
    end

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC}, eval::AbstractEvaluator, params::AbstractDict)
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N / T * (mag2 - mag^2)
    end
    evaluate!(eval, :BinderRatio, (:Mag2, :Mag4)) do mag2, mag4
        1 - (mag4/3/mag2^2)
    end

    C3 = [-1/2 -√3/2 0; √3/2 1/2 0; 0 0 1]
    for f in corr_posns
        evaluate!(eval, Symbol("χs_", f), (Symbol("sk_", f), Symbol("sk_corr_", f))) do sk, sk2
            N / T * (sk2 - sum(abs2.(sk)))
        end
        evaluate!(eval, Symbol("χη_", f), (Symbol("ηk_", f), Symbol("ηk_corr_", f))) do ηk, ηk2
            N / T * abs.(ηk2 - ηk*ηk')
        end
        if f in (M, part_K, half_M)
            corrnames = (Symbol("sk_corr_", f), Symbol("sk_corr_", f, "2"), Symbol("sk_corr_", f, "3"))
            evaluate!(eval, Symbol("sk_corr_", f, "_c3"), corrnames) do sk1, sk2, sk3
                sk1 + sk2 + sk3
            end
            corrnames = (Symbol("ηk_corr_", f), Symbol("ηk_corr_", f, "2"), Symbol("ηk_corr_", f, "3"))
            evaluate!(eval, Symbol("ηk_corr_", f, "_c3"), corrnames) do ηk1, ηk2, ηk3
                ηk1 + C3' * ηk2 * C3 + C3 * ηk3 * C3'
            end
        end
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end
function Carlo.register_evaluables(::Type{WignerMC{AlgType, BiasType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType, BiasType}
    Carlo.register_evaluables(WignerMC, eval, params)
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
