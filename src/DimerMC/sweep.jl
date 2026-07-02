function sweep_η!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng
    for _ in 1:length(mc.spins)
        # Select site for spin change
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

        old_E = site_energy(mc, mc.ηs[pos...], pos)
        new_η = rand(rng, SpinVector)
        new_E = site_energy(mc, new_η, pos)
        ΔE = new_E - old_E

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / mc.T)
        if prob >= 1.0 || rand(rng) < prob
            mc.ηs[pos...] = new_η
        end
    end
end

function sweep_s!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng

    changed = 0
    for _ in 1:max(Lx, Ly)
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    end

    if is_thermalized(ctx)
        measure!(ctx, :DimerChanges, changed)
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    if !mc.etaonly
        sweep_s!(mc, ctx)
    end
    sweep_η!(mc, ctx)
end