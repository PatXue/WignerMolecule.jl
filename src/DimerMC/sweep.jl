function metropolisacc(mc, ctx::Carlo.MCContext, ΔE)
    return ΔE <= 0.0 || rand(ctx.rng) < exp(-ΔE / calc_temp(mc, ctx))
end

function sweep_η!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng

    for _ in 1:length(mc.spins)
        # Select site for spin change
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

        old_E = site_energy_eta(mc, pos, mc.ηs[pos...])
        new_η = rand(rng, SpinVector)
        new_E = site_energy_eta(mc, pos, new_η)
        ΔE = new_E - old_E

        if metropolisacc(mc, ctx, ΔE)
            mc.ηs[pos...] = new_η
        end
    end
end

function sweep_dimer!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    rng = ctx.rng
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    posj = pos + rand(rng, disps)

    if mod_equiv(mc.spins[pos...], posj, mc) # Dimer dissolution
        s = rand(rng, SpinVector)
        sj = rand(rng, SpinVector)
        old_E = dimer_energy_s(mc, Dimer(pos, posj))
        new_E = pair_energy_s(mc, pos, posj, s, sj)
        if metropolisacc(mc, ctx, new_E - old_E)
            addmonomer!(pos, s, mc)
            addmonomer!(posj, sj, mc)
        end
    elseif ismonomer(pos, mc) && ismonomer(posj, mc) # Dimer creation
        s = mc.monospins[pos...]
        sj = mc.monospins[posj...]
        old_E = pair_energy_s(mc, pos, posj, s, sj)
        new_E = dimer_energy_s(mc, Dimer(pos, posj))
        if metropolisacc(mc, ctx, new_E - old_E)
            delmonomer!(pos, mc)
            delmonomer!(posj, mc)
            mc.spins[pos...] = posj
            mc.spins[posj...] = pos
        end
    elseif ismonomer(pos, mc) || ismonomer(posj, mc) # Monomer shift
        if ismonomer(posj, mc)
            pos, posj = posj, pos
        end
        s = mc.monospins[pos...]
        sk = rand(rng, SpinVector)
        posk = mc.spins[posj...]
        old_E = shift_energy_s(mc, Dimer(posj, posk), pos, s)
        new_E = shift_energy_s(mc, Dimer(pos, posj), posk, sk)
        if metropolisacc(mc, ctx, new_E - old_E)
            delmonomer!(pos, mc)
            addmonomer!(posk, sk, mc)
            mc.spins[pos...] = posj
            mc.spins[posj...] = pos
        end
    end
end

function sweep_monomer!(mc::DimerMC, ctx::Carlo.MCContext)
    rng = ctx.rng
    pos = randmonomer(mc, rng)

    old_E = site_energy_s(mc, pos, mc.spins[pos...])
    new_s = rand(rng, SpinVector)
    new_E = site_energy_s(mc, pos, new_s)
    ΔE = new_E - old_E

    if metropolisacc(mc, ctx, ΔE)
        mc.spins[pos...] = new_s
    end
end

function sweep_s!(mc::DimerMC, ctx::Carlo.MCContext)
    for _ in 1:length(mc.spins)
        if rand(ctx.rng) < mc.Q * length(mc.monomers) / length(mc.spins)
            sweep_monomer!(mc, ctx)
        else
            sweep_dimer!(mc, ctx)
        end
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    if !mc.etaonly
        sweep_s!(mc, ctx)
    end
    sweep_η!(mc, ctx)
end