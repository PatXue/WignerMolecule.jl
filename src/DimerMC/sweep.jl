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

function sweep_s!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    T = calc_temp(mc, ctx)
    rng = ctx.rng

    for _ in 1:length(mc.spins)
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
        posj = pos + rand(rng, disps)
        if mod_equiv(mc.spins[pos...], posj, mc)         # Dimer dissolution
            s = rand(SpinVector)
            sj = rand(SpinVector)
            old_E = dimer_energy_s(mc, Dimer(pos, posj))
            new_E = pair_energy_s(mc, pos, posj, s, sj)
        elseif ismonomer(pos, mc) && ismonomer(posj, mc) # Dimer creation
            s = mc.spins[pos...]
            sj = mc.spins[posj...]
            old_E = pair_energy_s(mc, pos, posj, s, sj)
            new_E = dimer_energy_s(mc, Dimer(pos, posj))
        elseif ismonomer(pos, mc) || ismonomer(posj, mc) # Monomer shift
            if ismonomer(posj, mc)
                pos, posj = posj, pos
            end
            s = rand(SpinVector)
            posk = mc.spins[posj...]
            old_E = shift_energy_s(mc, Dimer(posj, posk), pos, mc.spins[pos...])
            new_E = shift_energy_s(mc, Dimer(pos, posj), posk, s)
        else                                             # Adjacent dimers
            continue
        end
    end

    if is_thermalized(ctx)
        measure!(ctx, :dimersteps, steps)
        measure!(ctx, :dimerretries, retries)
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    if !mc.etaonly
        sweep_s!(mc, ctx)
    end
    sweep_η!(mc, ctx)
end