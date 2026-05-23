function Carlo.measure!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly

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
    for f in (Γ, M, M2, M3, half_M, half_K)
        η = f(mc.ηks)
        measure!(ctx, Symbol("ηk_", f), η)
        measure!(ctx, Symbol("ηk_corr_", f), η*η')
    end

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{DimerMC}, eval::AbstractEvaluator, params::AbstractDict)
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :ChiEtaxy, (:ηk_corr_Γ, :ηk_Γ)) do eta2, eta
        return N * real(eta2[1,1] + eta2[2,2] - abs2(eta[1]) - abs2(eta[2])) / T
    end
    evaluate!(eval, :ChiEtaz, (:ηk_corr_Γ, :ηk_Γ)) do eta2, eta
        return N * real(eta2[3,3] - abs2(eta[3])) / T
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end

function Carlo.write_checkpoint(mc::DimerMC, out::HDF5.Group)
    out["spins"] = Int.(mc.spins)
    out["etas"] = mc.ηs
    return nothing
end
function Carlo.read_checkpoint!(mc::DimerMC, in::HDF5.Group)
    mc.spins .= Bond.(read(in, "spins"))
    raw_ηs = read(in, "etas")
    mc.ηs .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_ηs)
    return nothing
end

