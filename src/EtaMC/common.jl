# Initialize the WignerMC spins
function Carlo.init!(mc::EtaMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = get(params, :init_type, :rand)

    rand!(ctx.rng, mc.spins)
    if init_type == :fm
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :fe
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 1, 0)
        end
    elseif init_type == :stripe
        for I in eachindex(mc.spins)
            _, y = Tuple(I)
            mc.spins[I] = SpinVector((-1)^y, 0, 0)
        end
    elseif init_type == :afm
        for I in eachindex(mc.spins)
            x, _ = Tuple(I)
            mc.spins[I] = SpinVector(0, 0, (-1)^x)
        end
    end

    update_fourier!(mc)
    return nothing
end

function Carlo.measure!(mc::EtaMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Site spin expectations
    mag_v = sum(mc.spins) ./ N
    mag = norm(mag_v)
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Magx, mag_v[1])
    measure!(ctx, :Magy, mag_v[2])
    measure!(ctx, :Magz, mag_v[3])

    # Energy per lattice site
    E = total_energy(mc) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in (Γ, M, M2, M3)
        pos = f(Lx, Ly)
        s = mc.spinks[pos..., :]
        χ = mc.chis[pos...]
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), s*s')
        measure!(ctx, Symbol("chik_", f), χ)
        measure!(ctx, Symbol("chik_corr_", f), abs2(χ))
    end

    single_Q = chirality(mc.spins, 1, 1)
    measure!(ctx, :single_Q, single_Q)
    measure!(ctx, :single_Q2, single_Q^2)

    return nothing
end

function Carlo.register_evaluables(::Type{EtaMC}, eval::AbstractEvaluator,
                                   params::AbstractDict)
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end
    return nothing
end

function Carlo.write_checkpoint(mc::EtaMC, out::HDF5.Group)
    out["spins"] = mc.spins
    return nothing
end
function Carlo.read_checkpoint!(mc::EtaMC, in::HDF5.Group)
    raw_spins = read(in, "spins")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_spins)
    return nothing
end
