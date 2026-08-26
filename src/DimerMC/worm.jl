function worm_dimer!(mc::DimerMC{:Worm}, ctx::Carlo.MCContext)
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

