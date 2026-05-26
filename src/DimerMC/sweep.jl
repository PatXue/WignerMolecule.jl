function sweep_η!(mc::DimerMC, ctx::Carlo.MCContext)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

        old_E = site_energy(mc, mc.ηs[pos...], pos)
        new_η = rand(rng, SpinVector)
        new_E = site_energy(mc, new_η, pos)
        ΔE = new_E - old_E

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / T)
        if prob >= 1.0 || rand(rng) < prob
            mc.ηs[x, y] = new_η
        end
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
end