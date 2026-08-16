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

    η = sum(mc.ηs) ./ N
    measure!(ctx, :ηz, abs(η[3]))
    measure!(ctx, :ηxy, sqrt(η[1]^2 + η[2]^2))

    # Energy per lattice site
    E = total_energy(mc) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in (Γ, M, half_M, part_K)
        pos = convert(SVector{2,Int}, f(Lx, Ly))
        s = getc3fourier(mc.spinks, pos)
        eta = getc3fourier(mc.ηks, pos)
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), norm2(s))
        measure!(ctx, Symbol("sk_quar_", f), norm2(s)^2)
        measure!(ctx, Symbol("etak_", f), eta)
        measure!(ctx, Symbol("etak_corr_", f), eta * eta')
        if mc.corr_rad != 0
            s = getc3fourier(mc.spinks, pos, mc.corr_rad)
            eta = getc3fourier(mc.ηks, pos, mc.corr_rad)
            measure!(ctx, Symbol("sk_near_", f), s)
            measure!(ctx, Symbol("sk_corr_near_", f), norm2(s))
            measure!(ctx, Symbol("sk_quar_near_", f), norm2(s)^2)
            measure!(ctx, Symbol("etak_near_", f), eta)
            measure!(ctx, Symbol("etak_corr_near_", f), eta * eta')
        end
    end

    for phase in (:fm, :stripe, :afm_fe, :afm_afe)
        if phase == :fm
            pos = SVector(1,1)
            a = SVector(0,0,1)
        elseif phase == :stripe
            pos = SVector{2,Int}(M(Lx, Ly))
            a = SVector(1/2,√3/2,0)
        elseif phase == :afm_fe
            pos = SVector(1,1)
            a = SVector(1/2,√3/2,0)
        elseif phase == :afm_afe
            pos = SVector{2,Int}(M2(Lx, Ly))
            a = SVector(0,1,0)
        end
        etak = a ⋅ getc3fourier(mc.ηks, pos)
        measure!(ctx, Symbol("etak_", phase), etak)
        measure!(ctx, Symbol("etak_corr_", phase), abs2(etak))
        measure!(ctx, Symbol("etak_quar_", phase), abs2(etak)^2)
        if mc.corr_rad != 0
            etak = a ⋅ getc3fourier(mc.ηks, pos, mc.corr_rad)
            measure!(ctx, Symbol("etak_near_", phase), etak)
            measure!(ctx, Symbol("etak_corr_near_", phase), abs2(etak))
            measure!(ctx, Symbol("etak_quar_near_", phase), abs2(etak)^2)
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

    for f in (Γ, M, half_M, part_K)
        evaluate!(eval, Symbol("χs_", f), (Symbol("sk_", f), Symbol("sk_corr_", f))) do sk, sk2
            N / T * (sk2 - sum(abs2.(sk)))
        end
        evaluate!(eval, Symbol("sk_kurt_", f), (Symbol("sk_corr_", f), Symbol("sk_quar_", f))) do sk2, sk4
            1 - sk4 / 3sk2^2
        end
        if get(params, :corr_rad, 0) != 0
            evaluate!(eval, Symbol("χs_near_", f), (Symbol("sk_near_", f), Symbol("sk_corr_near_", f))) do sk, sk2
                N / T * (sk2 - sum(abs2.(sk)))
            end
            evaluate!(eval, Symbol("sk_kurt_near_", f), (Symbol("sk_corr_near_", f), Symbol("sk_quar_near_", f))) do sk2, sk4
                1 - sk4 / 3sk2^2
            end
        end
    end

    for phase in (:fm, :stripe, :afm_fe, :afm_afe)
        evaluate!(eval, Symbol("χeta_", phase), (Symbol("etak_", phase), Symbol("etak_corr_", phase))) do sk, sk2
            N / T * (sk2 - sum(abs2.(sk)))
        end
        evaluate!(eval, Symbol("etak_kurt_", phase), (Symbol("etak_corr_", phase), Symbol("etak_quar_", phase))) do sk2, sk4
            1 - sk4 / 3sk2^2
        end
        if get(params, :corr_rad, 0) != 0
            evaluate!(eval, Symbol("χeta_near_", phase), (Symbol("etak_near_", phase), Symbol("etak_corr_near_", phase))) do sk, sk2
                N / T * (sk2 - sum(abs2.(sk)))
            end
            evaluate!(eval, Symbol("etak_kurt_near_", phase), (Symbol("etak_corr_near_", phase), Symbol("etak_quar_near_", phase))) do sk2, sk4
                1 - sk4 / 3sk2^2
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
