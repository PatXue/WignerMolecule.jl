function Carlo.sweep!(mc::MC{:Metropolis}, rng::AbstractRNG=default_rng())
    Lx, Ly = size(mc.spins)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        x = rand(rng, 1:Lx)
        y = rand(rng, 1:Ly)
        old_s = mc.spins[x, y]
        # Propose new spin vector
        new_s = rand(rng, SpinVector)
        ΔE = energy(mc, new_s, x, y) - energy(mc, old_s, x, y)

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / mc.T)
        if prob >= 1.0 || rand(rng) < prob
            mc.spins[x, y] = new_s
        end
    end
    return nothing
end
