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
    fill!(mc.visited, false)

    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    offset = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    rotation = rand(rng, 0:5)
    reflect_type = rand(rng, Bool)

    old_E = new_E = 0.0
    pocket::Vector{Dimer} = [Dimer(pos, mc.spins[pos...])]
    proposal::Vector{Dimer} = [] # Proposed new dimers
    while length(pocket) > 0
        d = pop!(pocket)
        mc.visited[d.pos...] = true
        mc.visited[d.posj...] = true
        old_E += bond_energy(mc, d)

        d = shift(d, -offset)
        d = invrotate(d, rotation)
        d = reflect_type ? reflect1(d) : reflect2(d)
        d = rotate(d, rotation)
        d = shift(d, offset)
        d = mod1(d, mc)
        push!(proposal, d)
        new_E += bond_energy(mc, d)
        append!(pocket, collisions(mc, d))
    end

    ΔE = new_E - old_E
    if ΔE < 0 || rand(rng) < exp(-ΔE / mc.T)
        for d in proposal
            mc.spins[d.pos...] = d.posj
            mc.spins[d.posj...] = d.pos
        end
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    sweep_s!(mc, ctx)
    sweep_η!(mc, ctx)
end