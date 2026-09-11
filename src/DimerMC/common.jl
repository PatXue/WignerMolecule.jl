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
        pos = convert(SVector{2,Int}, f(Lx, Ly))
        if mc.corr_rad == 0
            s = mc.sks[pos..., :]
            scorr = s * s'
            eta = mc.ηks[pos..., :]
            etacorr = eta * eta'
        else
            x, y = pos[1], pos[2]
            r = mc.corr_rad
            s = sum(eachslice(mc.sks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
            scorr = sum(sk -> sk * sk', mc.sks[x-r:x+r, y-r:y+r, :])
            eta = sum(eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
            etacorr = sum(etak -> etak * etak', eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
        end
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), scorr)
        measure!(ctx, Symbol("etak_", f), eta)
        measure!(ctx, Symbol("etak_corr_", f), etacorr)
    end

    for phase in (:fm, :stripe, :afm_fe, :afm_afe)
        if phase == :fm
            posns = [SVector(1,1)]
            as = [SVector(0.0,0,1)]
        elseif phase == :stripe
            posns = [SVector{2,Int}(M(Lx, Ly)), SVector{2,Int}(M2(Lx, Ly)), SVector{2,Int}(M3(Lx, Ly))]
            as = [SVector(1/2,√3/2,0), SVector(-1.0,0,0), SVector(1/2,-√3/2,0)]
        elseif phase == :afm_fe
            posns = [SVector(1,1),SVector(1,1),SVector(1,1)]
            as = [SVector(1/2,√3/2,0), SVector(-1.0,0,0), SVector(1/2,-√3/2,0)]
        elseif phase == :afm_afe
            posns = [SVector{2,Int}(M2(Lx, Ly)), SVector{2,Int}(M3(Lx, Ly)), SVector{2,Int}(M(Lx, Ly))]
            as = [SVector(0.0,1,0), SVector(-√3/2,-1/2,0), SVector(√3/2,-1/2,0)]
        end
        etatot = 0.0 + 0.0im
        etacorr = 0.0
        for (pos, a) in Iterators.zip(posns, as)
            if mc.corr_rad != 0
                x, y = pos[1], pos[2]
                r = mc.corr_rad
                etatot += sum(etak -> a ⋅ etak, eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
                etacorr += sum(etak -> abs2(a ⋅ etak), eachslice(mc.ηks[x-r:x+r, y-r:y+r, :], dims=(1,2)))
            else
                etak = mc.ηks[pos..., :]
                etatot += a ⋅ etak
                etacorr += abs2(a ⋅ etak)
            end
        end
        measure!(ctx, Symbol("etak_re_", phase), real(etatot))
        measure!(ctx, Symbol("etak_im_", phase), imag(etatot))
        measure!(ctx, Symbol("etak_corr_", phase), etacorr)
        measure!(ctx, Symbol("etak_quar_", phase), etacorr^2)
    end

    mc.sks .= abs2.(mc.sks)
    mc.ηks .= abs2.(mc.ηks)
    ifft!(mc.sks.array, (1,2))
    ifft!(mc.ηks.array, (1,2))
    for j in 1:3
        sr_corrs = zeros(div(Lx,2), 4)
        ηr_corrs = zeros(div(Lx,2), 3)
        a = oriented_disps[j]
        for i in 0:(div(Lx,2)-1)
            sr_corrs[i+1,:] .+= real.(mc.sks[mod1.([1,1] + i*a, (Lx, Ly))..., :])
            ηr_corrs[i+1,:] .+= real.(mc.ηks[mod1.([1,1] + i*a, (Lx, Ly))..., :])
            sr_corrs[i+1,:] .+= real.(mc.sks[mod1.([1,1] - i*a, (Lx, Ly))..., :])
            ηr_corrs[i+1,:] .+= real.(mc.ηks[mod1.([1,1] - i*a, (Lx, Ly))..., :])
        end
        sr_corrs ./= 2
        ηr_corrs ./= 2
        measure!(ctx, Symbol("sr_corr_a$j"), sr_corrs)
        measure!(ctx, Symbol("etar_corr_a$j"), ηr_corrs)
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

