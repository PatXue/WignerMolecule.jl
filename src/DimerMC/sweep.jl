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
        prob = exp(-ΔE / T)
        if prob >= 1.0 || rand(rng) < prob
            mc.ηs[x, y] = new_η
        end
    end
end

function sweep_s!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng

    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    offset = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    rotation = Bond(rand(rng, 0:5))
    reflect_type = rand(rng, Bool)

    mc.spins[pos...] = none
    mc.spins[getpartner(mc, pos)...] = none

    pocket::Vec{Dimer} = [Dimer(pos, getpartner(mc, pos))]
    while length(pocket) > 0
        d = pop!(pocket)
        d = shift(d, -offset)
        d = invrotate(d, rotation)
        d = reflect_type ? reflect1(d) : reflect2(d)
        d = rotate(d, rotation)
        d = shift(d, offset)
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
end