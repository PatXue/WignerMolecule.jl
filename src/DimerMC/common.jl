function Carlo.measure!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly

    η = sum(mc.ηs) ./ N
    measure!(ctx, :ηx, η[1])
    measure!(ctx, :ηy, η[2])
    measure!(ctx, :ηz, abs(η[3]))
    measure!(ctx, :ηxy, sqrt(η[1]^2 + η[2]^2))

    # Energy per lattice site
    E = total_energy(mc) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in corr_posns
        pos = f(Lx,Ly)
        s = mc.sks[pos...]
        η = mc.ηks[pos..., :]
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), abs2(s))
        measure!(ctx, Symbol("ηk_", f), η)
        measure!(ctx, Symbol("ηk_corr_", f), η*η')
    end
    return nothing
end

function Carlo.register_evaluables(::Type{DimerMC}, eval::AbstractEvaluator, params::AbstractDict)
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end
    return nothing
end

function Carlo.write_checkpoint(mc::DimerMC, out::HDF5.Group)
    out["spins"] = mc.spins
    out["spinks"] = mc.spinks
    out["etas"] = mc.ηs
    out["etaks"] = mc.spinks
    return nothing
end
function Carlo.read_checkpoint!(mc::DimerMC, in::HDF5.Group)
    raw_ss = read(in, "spins")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2]), raw_ss)
    raw_ηs = read(in, "etas")
    mc.ηs .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_ηs)
    return nothing
end

