function sweep_η!(mc::DimerMC, T, rng=default_rng())
    Lx, Ly = size(mc.ηs)
    for _ in 1:length(mc.ηs)
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
        mc.ηs[pos...] = heatbath_spin(site_field_eta(mc, pos), T; rng)
    end
end

function sweep_dimer!(mc::DimerMC, T, rng=default_rng())
    Lx, Ly = size(mc.spins)

    for _ in 1:length(mc.spins)
        pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
        posj = pos + rand(rng, disps)

        if mod_equiv(mc.spins[pos...], posj, mc) # Dimer dissolution
            s = rand(rng, SpinVector)
            sj = rand(rng, SpinVector)
            old_E = dimer_energy_s(mc, Dimer(pos, posj))
            new_E = pair_energy_s(mc, pos, posj, s, sj)
            if metropolisacc(new_E - old_E, T; rng, fug=1/mc.fug)
                addmonomer!(pos, s, mc)
                addmonomer!(posj, sj, mc)
            end

        elseif ismonomer(pos, mc) && ismonomer(posj, mc) # Dimer creation
            s = mc.monospins[pos...]
            sj = mc.monospins[posj...]
            old_E = pair_energy_s(mc, pos, posj, s, sj)
            new_E = dimer_energy_s(mc, Dimer(pos, posj))
            if metropolisacc(new_E - old_E, T; rng, fug=mc.fug)
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
            if metropolisacc(new_E - old_E, T; rng)
                delmonomer!(pos, mc)
                addmonomer!(posk, sk, mc)
                mc.spins[pos...] = posj
                mc.spins[posj...] = pos
            end

        else # Double dimer move
            d = orientdimer(getdimer(pos, mc), mc)
            dj = orientdimer(getdimer(posj, mc), mc)
            if arestacked(d, dj, mc)
                new_d = Dimer(d.pos, dj.pos)
                new_dj = Dimer(d.posj, dj.posj)
                old_E = dimer_energy_s(mc, d) + dimer_energy_s(mc, dj)
                new_E = dimer_energy_s(mc, new_d) + dimer_energy_s(mc, new_dj)
                if metropolisacc(new_E - old_E, T; rng)
                    add_dimer!(new_d, mc)
                    add_dimer!(new_dj, mc)
                end
            end
        end
    end
end

function worm_dimer!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    T = calc_temp(mc, ctx)
    copy!(mc.spinscopy, mc.spins)
    rng = ctx.rng

    steps = 0
    retries = 0
    final_pos = pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    while true
        Zs = [exp(-bond_energy(mc, Dimer(pos, pos+a)) / T) for a in disps]
        Z = sum(Zs)
        p = rand(rng)
        for i in 1:6
            p -= Zs[i] / Z
            if i == 6 || p < 0
                posj = pos + disps[i]
                if !mod_equiv(mc.spins[pos...], posj, mc)
                    steps += 1
                end
                mc.spins[pos...] = posj
                pos, mc.spins[posj...] = mc.spins[posj...], pos
                break
            end
        end
        if mod_equiv(final_pos, pos, mc)
            break
        elseif steps > 100 * Lx * Ly    # Retry after loop exceeds 100N
            steps = 0
            retries += 1
            final_pos = pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
            copy!(mc.spins, mc.spinscopy)
        end
    end
end

function sweep_monomer!(mc::DimerMC, T, rng=default_rng())
    for _ in 1:length(mc.monomers)
        pos = randmonomer(mc, rng)

        old_E = site_energy_s(mc, pos, mc.monospins[pos...])
        new_s = rand(rng, SpinVector)
        new_E = site_energy_s(mc, pos, new_s)

        if metropolisacc(new_E - old_E, T; rng)
            mc.monospins[pos...] = new_s
        end
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    sweep_dimer!(mc, calc_temp(mc, ctx), ctx.rng)
    sweep_monomer!(mc, calc_temp(mc, ctx), ctx.rng)
    sweep_η!(mc, calc_temp(mc, ctx), ctx.rng)
end

function Carlo.sweep!(mc::DimerMC{:SpinOnly}, ctx::Carlo.MCContext)
    sweep_dimer!(mc, calc_temp(mc, ctx), ctx.rng)
    sweep_monomer!(mc, calc_temp(mc, ctx), ctx.rng)
end

function Carlo.sweep!(mc::DimerMC{:HeatbathEta}, ctx::Carlo.MCContext)
    sweep_η!(mc, calc_temp(mc, ctx), ctx.rng)
end
