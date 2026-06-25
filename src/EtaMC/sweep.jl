function Carlo.sweep!(mc::EtaMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng
    T = calc_temp(mc, ctx)
    for _ in 1:length(mc.spins)
        # Select site for spin change
        x = rand(rng, 1:Lx)
        y = rand(rng, 1:Ly)

        old_s = mc.spins[x, y]
        new_s = rand(rng, SpinVector)
        ΔE = energy(mc, new_s - old_s, x, y)    # energy linear in s

        # Probability of accepting spin flip (for ΔE ≤ 0 always accept)
        prob = exp(-ΔE / T)
        if prob >= 1.0 || rand(rng) < prob
            mc.spins[x, y] = new_s
        end
    end
    return nothing
end
