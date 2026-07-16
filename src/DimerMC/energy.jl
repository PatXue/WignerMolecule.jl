#=
Helper functions for calculating system energy
Broken into sections for the η sweeps, s sweeps, and total energy calculations

Functions for η sweeps take ηs as arguments and use the mc's current spin state
Functions for s sweeps take the spin state as arguments and use mc's η state,
    they also only calculate the spin-orbit coupling energy
Functions for total energy are used for measurements, so only use the current mc state
=#

## η sweep functions ##

"""
    ssfactor(mc::DimerMC, η, ηj, ν)

Calculate `s⋅sj` coefficient given `(η, ηj, ν)` in Hamiltonian. `η`s expected to be half-unit vectors.
"""
function ssfactor(mc::DimerMC, η, ηj, ν)
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS

    # η raising and lowering operators
    η_m = η[1] + 1.0im*η[2]
    ηj_p = ηj[1] - 1.0im*ηj[2]
    ηj_m = ηj[1] + 1.0im*ηj[2]

    E_spin = 0.0 + 0.0im
    E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
    E_spin += 2*J_EMEP_SS *     η_m * ηj_p
    E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
    E_spin += J_SS
    E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)
    return real(E_spin)
end

function get_sdot(d::Dimer, mc::DimerMC)
    if mod_equiv(mc.spins[d.pos...], d.posj, mc)
        return -3/4
    elseif ismonomer(d.pos, mc) && ismonomer(d.posj, mc)
        return mc.monospins[d.pos...] ⋅ mc.monospins[d.posj...] / 4
    else
        return 0.0
    end
end

function bond_energy(mc::DimerMC, d::Dimer, η, ηj)
    # Couplings
    J_EzEz = mc.params.J_EzEz
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    ν = getν(d, mc)
    sdot = get_sdot(d, mc)
    # η raising and lowering operators
    η_m = η[1] + 1.0im*η[2]
    ηj_p = ηj[1] - 1.0im*ηj[2]
    ηj_m = ηj[1] + 1.0im*ηj[2]

    E_spin = 0.0
    E_η = 0.0 + 0.0im

    # η-only energy
    E_η +=   J_EzEz *     η[3] * ηj[3]
    E_η += 2*J_EMEP *     η_m * ηj_p
    E_η += 2*J_EMEM * ν * η_m * ηj_m

    # η-S energy
    if sdot != 0.0
        E_spin = sdot * ssfactor(mc, η, ηj, ν)
    end

    return E_spin + real(E_η)
end

function site_energy_eta(mc::DimerMC, pos, η)
    η /= 2
    E = 0.0
    for j in 1:3
        disp = oriented_disps[j]
        posj = pos + disp
        ηj = mc.ηs[posj...] / 2
        E += bond_energy(mc, Dimer(pos, posj), η, ηj)
        posj = pos - disp
        ηj = mc.ηs[posj...] / 2
        E += bond_energy(mc, Dimer(posj, pos), ηj, η)
    end
    return E
end

## Spin sweep functions ##

"""
    ssfactor(mc::DimerMC, d::Dimer)

Calculate `s⋅sj` coefficient for a bond given by `d`
"""
function ssfactor(mc::DimerMC, d::Dimer)
    d = orientdimer(d, mc)
    ν = getν(d, mc)
    η = mc.ηs[d.pos...] / 2
    ηj = mc.ηs[d.posj...] / 2
    return ssfactor(mc, η, ηj, ν)
end

# Energy from spin-orbit coupling on bond d with given sdot
bond_energy_s(mc::DimerMC, d::Dimer, sdot) = sdot * ssfactor(mc, d)

# Spin orbit energy of a single monomer and entangled dimer (when shifting a monomer)
function shift_energy_s(mc::DimerMC, d::Dimer, pos, s)
    E = bond_energy_s(mc, d, -3/4)
    for disp in disps
        posj = pos + disp
        if !ismonomer(posj, mc) || indimer(pos, d, mc)
            continue
        else
            sdot = s ⋅ mc.monospins[posj...]
        end
        E += bond_energy_s(mc, Dimer(pos, posj), sdot)
    end
    return E
end

# Spin-orbit energy of a pair of adjacent monomers (when dissolving/forming a dimer)
function pair_energy_s(mc::DimerMC, d::Dimer, s, sj)
end

## Total energy functions ##

function bond_energy(mc::DimerMC, d::Dimer)
    d = orientdimer(d, mc)
    return bond_energy(mc, d, mc.ηs[d.pos...], mc.ηs[d.posj...])
end

# Energy from half the bonds of pos
function half_energy(mc::DimerMC, pos)
    E = 0.0
    for disp in oriented_disps
        posj = pos .+ disp
        E += bond_energy(mc, Dimer(pos, posj))
    end
    return E
end

function total_energy(mc::DimerMC)
    E = 0.0
    for I in eachindex(IndexCartesian(), mc.spins)
        E += half_energy(mc, Tuple(I))
    end
    return E
end

