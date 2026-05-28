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
    E = total_energy(mc, mc.B) / N
    measure!(ctx, :Energy, E)
    measure!(ctx, :Energy2, E^2)

    update_fourier!(mc)
    for f in (Γ, M, M2, M3)
        s = f(mc.spinks)
        measure!(ctx, Symbol("sk_", f), s)
        measure!(ctx, Symbol("sk_corr_", f), s*s')
    end

    Q = chirality(mc.spins) / N
    measure!(ctx, :Q, Q)
    measure!(ctx, :Q2, Q^2)

    return nothing
end

