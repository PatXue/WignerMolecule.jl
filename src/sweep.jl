function metropolisacc(ΔE, T; rng, fug=1.0)
    prob = fug * exp(-ΔE/T)
    return prob >= 1.0 || rand(rng) <= prob
end

function flip_eta!(mc::WignerMC{:Metropolis}, T, rng)
    Lx, Ly = size(mc.ηs)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

    old_E = site_energy_eta(mc, pos, mc.ηs[pos...])
    new_η = rand(rng, SpinVector)
    new_E = site_energy_eta(mc, pos, new_η)
    ΔE = new_E - old_E

    if metropolisacc(ΔE, T; rng)
        mc.ηs[pos...] = new_η
    end
end

function flip_spin!(mc::WignerMC{:Metropolis}, T, B, rng)
    Lx, Ly = size(mc.ηs)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

    old_E = site_energy_s(mc, pos, mc.spins[pos...], B)
    new_s = rand(rng, SpinVector)
    new_E = site_energy_s(mc, pos, new_s, B)
    ΔE = new_E - old_E

    if metropolisacc(ΔE, T; rng)
        mc.spins[pos...] = new_s
    end
end

function Carlo.sweep!(mc, ctx::Carlo.MCContext)
    rng = ctx.rng
    T = calc_temp(mc, ctx)
    for _ in 1:length(mc.spins)
        flip_eta!(mc, T, rng)
        flip_spin!(mc, T, calc_B(mc, ctx), rng)
    end
    return nothing
end
