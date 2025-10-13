function Carlo.sweep!(mc::WignerMC{:Metropolis}, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng
    T = calc_therm_temp(mc, ctx)
    for _ in 1:length(mc.spins)
        for type in (:s, :η)
            # Select site for spin change
            x = rand(rng, 1:Lx)
            y = rand(rng, 1:Ly)

            old_s = mc.spins[x, y]
            old_η = mc.ηs[x, y]
            old_E = energy(mc, old_s, old_η, x, y)
            new_s = new_η = rand(rng, SpinVector)
            if type == :s
                new_E = energy(mc, new_s, old_η, x, y)
            elseif type == :η
                new_E = energy(mc, old_s, new_η, x, y)
            end
            ΔE = new_E - old_E

            # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
            prob = exp(-ΔE / T)
            if prob >= 1.0 || rand(rng) < prob
                if type == :s
                    mc.spins[x, y] = new_s
                elseif type == :η
                    mc.ηs[x, y] = new_η
                end
            end
        end
    end
    return nothing
end
