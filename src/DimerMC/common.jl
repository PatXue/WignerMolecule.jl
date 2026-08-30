function Carlo.measure!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly

    η = sum(mc.ηs) ./ N
    measure!(ctx, :ηz, abs(η[3]))
    measure!(ctx, :ηxy, sqrt(η[1]^2 + η[2]^2))

    # Energy per lattice site
    E = total_energy(mc) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in corr_posns
        pos = f(Lx,Ly)
        s = mc.sks[pos..., :]
        η = mc.ηks[pos..., :]
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), abs2.(s))
        measure!(ctx, Symbol("ηk_", f), η)
        measure!(ctx, Symbol("ηk_corr_", f), η*η')
        if mc.corr_rad != 0
            r = mc.corr_rad
            x, y = pos[1], pos[2]
            scorr = sum(sk -> abs2.(sk), mc.sks[x-r:x+r, y-r:y+r, :])
            ηcorr = sum(ηk -> ηk*ηk', eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
            measure!(ctx, Symbol("sk_corr_near_", f), scorr)
            measure!(ctx, Symbol("ηk_corr_near_", f), ηcorr)
        end
    end

    mc.sks .= abs2.(mc.sks)
    mc.ηks .= abs2.(mc.ηks)
    ifft!(mc.sks.array, (1,2))
    ifft!(mc.ηks.array, (1,2))
    sr_corrs = zeros(div(Lx,2), 4)
    ηr_corrs = zeros(div(Lx,2), 3)
    for a in disps
        for i in 0:(div(Lx,2)-1)
            sr_corrs[i+1,:] .+= real.(mc.sks[mod1.([1,1] + i*a, (Lx, Ly))..., :])
            ηr_corrs[i+1,:] .+= real.(mc.ηks[mod1.([1,1] + i*a, (Lx, Ly))..., :])
        end
    end
    sr_corrs ./= 6
    ηr_corrs ./= 6
    measure!(ctx, :sr_corrs, sr_corrs)
    measure!(ctx, :ηr_corrs, ηr_corrs)

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
    out["monospins"] = mc.monospins
    out["etas"] = mc.ηs
    out["Nw_val"] = (mc.Nw[]).val
    out["Nw_n"] = (mc.Nw[]).n
    return nothing
end
function Carlo.read_checkpoint!(mc::DimerMC, in::HDF5.Group)
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2]), read(in, "spins"))
    mc.monospins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), read(in, "monospins"))
    for I in eachindex(mc.spins, mc.monospins)
        pos = convert(SVector, I)
        if ismonomer(pos, mc)
            addmonomer!(pos, mc.monospins[I], mc)
        end
    end
    raw_ηs = read(in, "etas")
    mc.ηs .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_ηs)
    mc.Nw[] = Expectation(read(in, "Nw_val"), read(in, "Nw_n"))
    return nothing
end

