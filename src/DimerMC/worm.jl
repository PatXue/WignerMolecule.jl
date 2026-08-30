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
function worm_dimer!(mc::DimerMC{:Worm}, init_pos, T, rng)
    steps = 0
    changes = 0
    retries = 0

    pos = init_pos # Current worm head position
    posf = mc.spins[pos...] # Partner of init_pos, termination point
    # Store old config along worm path
    worm = Dict{Int, Int}(postoint(pos, mc) => postoint(posf, mc), postoint(posf, mc) => postoint(pos, mc))
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
        if ismonomer(posj, mc)              # Terminate at monomer
            delmonomer!(posj, mc)
            addmonomer!(posf, rand(rng, SpinVector), mc)
            mc.spins[posj...] = pos
            break
        elseif mod_equiv(posj, posf, mc)    # Terminate as loop
            mc.spins[posj...] = pos
            break
        end
        if !haskey(worm, postoint(posj, mc))
            nj = postoint(posj, mc)
            nk = postoint(mc.spins[posj...], mc)
            worm[nj] = nk
            worm[nk] = nj
        end
        pos, mc.spins[posj...] = mc.spins[posj...], pos
        steps += 1

        if steps > 100 * length(mc.spins) # Retry after loop exceeds 100N
            steps = 0
            retries += 1
            pos = randdimer(mc, rng)
            posf = mc.spins[pos...]
            for (n, nj) in pairs(worm)
                mc.spins[inttopos(n, mc)...] = inttopos(nj, mc)
            end
        end
    end
    return (changes, retries)
end

"""
    worm_mono!(mc::DimerMC{:Worm}, init_pos, ctx::Carlo.MCContext)

Compute a worm/loop update starting at the monomer site `init_pos`. Can
terminate at any point, forming a monomer at the head
"""
function worm_mono!(mc::DimerMC{:Worm}, init_pos, T, rng)
    steps = 0
    changes = 0
    retries = 0

    pos = init_pos # Current worm head position
    new_s = rand(rng, SpinVector)
    # Store old config along worm path
    worm = Dict{Int, Int}(postoint(pos, mc) => postoint(pos, mc))
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
        if !haskey(worm, postoint(posj, mc))
            nj = postoint(posj, mc)
            nk = postoint(mc.spins[posj...], mc)
            worm[nj] = nk
            worm[nk] = nj
        end
        mc.spins[pos...] = posj
        pos, mc.spins[posj...] = mc.spins[posj...], pos
        steps += 1

        if steps > 100 * length(mc.spins) # Retry after loop exceeds 100N
            steps = 0
            retries += 1
            pos = randmonomer(mc, rng)
            for (n, nj) in pairs(worm)
                mc.spins[inttopos(n, mc)...] = inttopos(nj, mc)
            end
        end
    end
    return (changes, retries)
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
                changes = worm_mono!(mc, pos, T, rng)[1]
            else
                changes = worm_dimer!(mc, pos, T, rng)[1]
            end
            tot_changes += changes + 1
            mc.Nw[] = addsample(mc.Nw[], changes+1)
        end
    else
        tot_changes = 0
        tot_retries = 0
        worms = cld(N, (mc.Nw[]).val)
        for _ in 1:worms
            pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
            if ismonomer(pos, mc)
                changes, retries = worm_mono!(mc, pos, T, rng)
            else
                changes, retries = worm_dimer!(mc, pos, T, rng)
            end
            tot_changes += changes
            tot_retries += retries
        end
        measure!(ctx, :WormChanges, tot_changes / worms)
        measure!(ctx, :WormRetries, tot_retries / worms)
    end
    sweep_monomer!(mc, T, ctx.rng)
    sweep_η!(mc, T, rng)
end
