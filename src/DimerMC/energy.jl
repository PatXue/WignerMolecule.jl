#=
Helper functions for calculating system energy
Broken into sections for the η and s sweeps calculations

Functions for η sweeps take ηs as arguments and use the mc's current spin state
Functions for s sweeps take the spin state as arguments and use mc's η state,
    they also only calculate the spin-orbit coupling energy
=#

## η sweep functions ##

function get_sdot(d::Dimer, mc::DimerMC)
    if mod_equiv(mc.spins[d.pos...], d.posj, mc)
        return -3/4
    elseif ismonomer(d.pos, mc) && ismonomer(d.posj, mc)
        return mc.monospins[d.pos...] ⋅ mc.monospins[d.posj...] / 4
    else
        return 0.0
    end
end

## Spin sweep functions ##

# Energy from spin-orbit coupling on bond d with given sdot
dimer_energy_s(mc::DimerMC, d::Dimer) = bond_energy_s(mc, d, -3/4)

function site_energy_s(mc::DimerMC, pos, s)
    E = 0.0
    for disp in disps
        posj = pos + disp
        if ismonomer(posj, mc)
            sdot = s ⋅ mc.monospins[posj...] / 4
            E += bond_energy_s(mc, Dimer(pos, posj), sdot)
        end
    end
    return E
end

# Spin orbit energy of a single monomer and entangled dimer (when shifting a monomer)
function shift_energy_s(mc::DimerMC, d::Dimer, pos, s)
    E = dimer_energy_s(mc, d)
    for disp in disps
        posj = pos + disp
        if ismonomer(posj, mc) && !indimer(posj, d, mc)
            sdot = s ⋅ mc.monospins[posj...] / 4
            E += bond_energy_s(mc, Dimer(pos, posj), sdot)
        end
    end
    return E
end

# Spin-orbit energy of a pair of adjacent monomers (when dissolving/forming a dimer)
function pair_energy_s(mc::DimerMC, pos, posj, s, sj)
    E = bond_energy_s(mc, Dimer(pos, posj), (s⋅sj) / 4)
    for disp in disps
        posk = pos + disp
        if !mod_equiv(posk, posj, mc) && ismonomer(posk, mc)
            sk = mc.monospins[posk...]
            E += bond_energy_s(mc, Dimer(pos, posk), (s⋅sk) / 4)
        end
        posk = posj + disp
        if !mod_equiv(posk, pos, mc) && ismonomer(posk, mc)
            sk = mc.monospins[posk...]
            E += bond_energy_s(mc, Dimer(posj, posk), (sj⋅sk) / 4)
        end
    end
    return E
end
