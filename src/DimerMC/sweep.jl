function sweep_η!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    T = calc_temp(mc, ctx)
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
            mc.ηs[pos...] = new_η
        end
    end
end

function sweep_s!(mc::DimerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    T = calc_temp(mc, ctx)
    rng = ctx.rng

    steps = 0
    changes = 0
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
                    changes += 1
                end
                mc.spins[pos...] = posj
                pos, mc.spins[posj...] = mc.spins[posj...], pos
                break
            end
        end
        steps += 1
        if mod_equiv(final_pos, pos, mc)
            break
        end
    end

    if is_thermalized(ctx)
        measure!(ctx, :dimersteps, steps)
        measure!(ctx, :dimerchanges, changes)
    end
end

function Carlo.sweep!(mc::DimerMC, ctx::Carlo.MCContext)
    if !mc.etaonly
        sweep_s!(mc, ctx)
        if ctx.sweeps % 1000 == 0
            println("Dimer sweep #$(ctx.sweeps) complete")
        end
    end
    sweep_η!(mc, ctx)
end