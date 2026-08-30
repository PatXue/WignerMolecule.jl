function sample(rng, as, Zs)
    Z = sum(Zs)
    prob = rand(rng)
    for i in eachindex(Zs, as)
        prob -= Zs[i] / Z
        if prob <= 0
            return as[i]
        end
    end
    return last(as)
end

"""
    worm_dimer!(mc::DimerMC{:Worm}, ctx::Carlo.MCContext)

Compute a worm/loop update starting at the dimer'd site `init_pos`. Can
end by forming a full loop, or by meeting a monomer site.
"""
function worm_dimer!(mc::DimerMC{:Worm}, init_pos, ctx::Carlo.MCContext)
    T = calc_temp(mc, ctx)
    copy!(mc.spinscopy, mc.spins)
    rng = ctx.rng

    steps = 0
    changes = 0
    retries = 0
    pos = init_pos # Current worm head position
    posf = mc.spins[pos...] # Partner of init_pos, termination point
    while true
        Zs = zeros(6)
        for i in 1:6
            posj = pos + disps[i]
            d = Dimer(pos, posj)
            if !ismonomer(posj, mc)
                Zs[i] = exp(-dimer_energy_s(mc, d) / T)
            else
                ΔE = dimer_energy_s(mc, d) - site_energy_s(mc, posj, mc.monospins[posj...])
                Zs[i] = exp(-ΔE / T)
            end
        end

        a = sample(rng, collect(disps), Zs)
        posj = pos + a
        if !mod_equiv(mc.spins[pos...], posj, mc)
            changes += 1
        end
        mc.spins[pos...] = posj
        if ismonomer(posj, mc)
            delmonomer!(posj, mc)
            addmonomer!(posf, rand(rng, SpinVector), mc)
            mc.spins[posj...] = pos
            break
        elseif mod_equiv(posj, posf, mc)
            mc.spins[posj...] = pos
            break
        end
        pos, mc.spins[posj...] = mc.spins[posj...], pos
        steps += 1

        if steps > 100 * length(mc.spins) # Retry after loop exceeds 100N
            steps = 0
            retries += 1
            pos = init_pos
            copy!(mc.spins, mc.spinscopy)
        end
    end
    return changes
end

"""
    worm_mono!(mc::DimerMC{:Worm}, init_pos, ctx::Carlo.MCContext)

Compute a worm/loop update starting at the monomer site `init_pos`. Can
terminate at any point, forming a monomer at the head
"""
function worm_mono!(mc::DimerMC{:Worm}, init_pos, ctx::Carlo.MCContext)
    T = calc_temp(mc, ctx)
    copy!(mc.spinscopy, mc.spins)
    rng = ctx.rng

    steps = 0
    changes = 0
    retries = 0
    pos = init_pos # Current worm head position
    new_s = rand(rng, SpinVector)
    while true
        Zs = zeros(7)
        for i in 1:6
            posj = pos + disps[i]
            if !ismonomer(posj, mc)
                d = Dimer(pos, posj)
                Zs[i] = exp(-dimer_energy_s(mc, d) / T)
            end
        end
        if sum(Zs) == 0.0
            break
        end
        Zs[7] = exp(-site_energy_s(mc, pos, new_s) / T)

        i = sample(rng, 1:7, Zs)
        if i == 7
            delmonomer!(init_pos, mc)
            addmonomer!(pos, new_s, mc)
            break
        end
        posj = pos + disps[i]
        if !mod_equiv(mc.spins[pos...], posj, mc)
            changes += 1
        end
        mc.spins[pos...] = posj
        pos, mc.spins[posj...] = mc.spins[posj...], pos
        steps += 1

        if steps > 100 * length(mc.spins) # Retry after loop exceeds 100N
            steps = 0
            retries += 1
            pos = init_pos
            copy!(mc.spins, mc.spinscopy)
        end
    end
    return changes
end

function Carlo.sweep!(mc::DimerMC{:Worm}, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    T = calc_temp(mc, ctx)
    rng = ctx.rng
    sweep_dimer!(mc, T, rng)
    if !is_thermalized(ctx)
        tot_changes = 0
        while tot_changes < N
            pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
            if ismonomer(pos, mc)
                changes = worm_mono!(mc, pos, ctx)
            else
                changes = worm_dimer!(mc, pos, ctx)
            end
            tot_changes += changes + 1
            mc.Nw[] = addsample(mc.Nw[], changes)
        end
    else
        tot_changes = 0
        for _ in 1:cld(N, (mc.Nw[]).val)
            pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
            if ismonomer(pos, mc)
                tot_changes += worm_mono!(mc, pos, ctx)
            else
                tot_changes += worm_dimer!(mc, pos, ctx)
            end
        end
        measure!(ctx, :WormChanges, tot_changes)
    end
    sweep_monomer!(mc, T, ctx.rng)
    sweep_η!(mc, T, rng)
end
