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
    addmonomer!(posf, rand(rng, SpinVector), mc)
    while true
        Zs = zeros(6)
        for i in 1:6
            posj = pos + disps[i]
            d = Dimer(pos, posj)
            if !ismonomer(posj, mc) || mod_equiv(posj, posf, mc)
                Zs[i] = exp(-dimer_energy_s(mc, d) / T)
            else
                Zs[i] = exp(-shift_energy_s(mc, d, posf, mc.monospins[posf...]) / T)
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
    return (steps, changes, retries)
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
    delmonomer!(pos, mc)
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
        Zs[7] = exp(-site_energy_s(mc, pos, new_s) / T)

        i = sample(rng, 1:7, Zs)
        if i == 7
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
    return (steps, changes, retries)
end

function Carlo.sweep!(mc::DimerMC{:Worm}, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    T = calc_temp(mc, ctx)
    rng = ctx.rng
    sweep_dimer!(mc, T, rng)
    for _ in 1:max(Lx, Ly)
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
        if ismonomer(pos, mc)
            worm_mono!(mc, pos, ctx)
        else
            worm_dimer!(mc, pos, ctx)
        end
    end
    sweep_monomer!(mc, T, ctx.rng)
end
