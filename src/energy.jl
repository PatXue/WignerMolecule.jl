#=
Helper functions for calculating system energy
Broken into sections for the η sweeps, s sweeps, and total energy calculations

Functions for η sweeps take ηs as arguments and use the mc's current spin state
Functions for s sweeps take the spin state as arguments and use mc's η state,
    they also only calculate the spin-orbit coupling energy
Functions for total energy are used for measurements, so only use the current mc state
=#

const ω::ComplexF64 = exp(im * 2π/3)

## η sweep functions ##

"""
    ssfactor(mc, η, ηj, ν)

Calculate `s⋅sj` coefficient given `(η, ηj, ν)` in Hamiltonian. `η`s expected to be half-unit vectors.
"""
function ssfactor(mc, η, ηj, ν)
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

get_sdot(d::Dimer, mc::WignerMC) = mc.spins[d.pos...] ⋅ mc.spins[d.posj...] / 4

function bond_energy(mc, d::Dimer, η, ηj)
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
    E_spin = sdot * ssfactor(mc, η, ηj, ν)

    return E_spin + real(E_η)
end

function site_energy_eta(mc, pos, η)
    η /= 2
    E = 0.0
    for disp in oriented_disps
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
    ssfactor(mc, d::Dimer)

Calculate `s⋅sj` coefficient for a bond given by `d`, which need not be oriented
"""
function ssfactor(mc, d::Dimer)
    d = orientdimer(d, mc)
    ν = getν(d, mc)
    η = mc.ηs[d.pos...] / 2
    ηj = mc.ηs[d.posj...] / 2
    return ssfactor(mc, η, ηj, ν)
end

# Energy from spin-orbit coupling on bond d with given sdot
bond_energy_s(mc, d::Dimer, sdot) = sdot * ssfactor(mc, d)

function site_energy_s(mc::WignerMC, pos, s)
    E = 0.0
    for disp in disps
        posj = pos + disp
        sdot = s ⋅ mc.spins[posj...] / 4
        E += bond_energy_s(mc, Dimer(pos, posj), sdot)
    end
    return E
end

function site_energy_s(mc::WignerMC{AlgType, Nothing}, pos, s, _) where AlgType
    return site_energy_s(mc, pos, s)
end
function site_energy_s(mc::WignerMC, pos, s, B)
    return site_energy_s(mc, pos, s) - B * mc.bias(pos...) ⋅ s
end

## Total energy functions ##

bond_energy(mc, d::Dimer) = bond_energy(mc, d, mc.ηs[d.pos...]/2, mc.ηs[d.posj...]/2)

# Energy from half the bonds of pos
function half_energy_nobias(mc, pos)
    E = 0.0
    for disp in oriented_disps
        posj = pos .+ disp
        E += bond_energy(mc, Dimer(pos, posj))
    end
    return E - mc.B * mc.bias(pos...) ⋅ mc.spins[pos...]
end
half_energy(mc, pos) = half_energy_nobias(mc, pos)
function half_energy(mc::WignerMC{AlgType, Nothing}, pos) where AlgType
    half_energy_nobias(mc, pos)
end
function half_energy(mc::WignerMC, pos)
    return half_energy_nobias(mc, pos) - mc.B * mc.bias(pos...) ⋅ mc.spins[pos...]
end

function total_energy(mc)
    E = 0.0
    for I in eachindex(IndexCartesian(), mc.spins)
        E += half_energy(mc, Tuple(I))
    end
    return E
end
