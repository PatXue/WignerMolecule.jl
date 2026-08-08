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

function bond_field_eta(mc, d::Dimer)
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EzEz = mc.params.J_EzEz
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    # For the site the field is calculated at (d.pos) use dummy vectors (as though η=(1,1,1)).
    # For misoriented d, swap η and ηj, so that energy formula is correct.
    # To avoid type instability, in both cases variables have the same type (may not be important)
    if isoriented(d, mc)
        η_z = SVector(0.0,0,1)
        η_m = SVector(1.0,im,0)
        η_p = SVector(1.0,-im,0)
        ηj = mc.ηs[d.posj...] / 2
        ηj_z = SVector(0.0,0,ηj[3])
        ηj_p = SVector(ηj[1] - 1.0im*ηj[2], ηj[1] - 1.0im*ηj[2], 0)
        ηj_m = conj.(ηj_p)
    else
        ηj_z = SVector(0.0,0,1)
        ηj_m = SVector(1.0,im,0)
        ηj_p = SVector(1.0,-im,0)
        η = mc.ηs[d.posj...] / 2
        η_z = SVector(0.0,0,η[3])
        η_p = SVector(η[1] - 1.0im*η[2], η[1] - 1.0im*η[2], 0)
        η_m = conj.(η_p)
    end
    ν = getν(orientdimer(d, mc), mc)
    sdot = get_sdot(d, mc)

    B = zeros(ComplexF64, 3)
    B_spin = zeros(ComplexF64, 3)

    # η-only energy
    B += J_EzEz   *     η_z .* ηj_z
    B += 2*J_EMEP *     η_m .* ηj_p
    B += 2*J_EMEM * ν * η_m .* ηj_m

    B_spin += J_EzEz_SS   *     η_z .* ηj_z
    B_spin += 2*J_EMEP_SS *     η_m .* ηj_p
    B_spin += 2*J_EMEM_SS * ν * η_m .* ηj_m
    B_spin += 2*J_EAM_SS  * (isoriented(d, mc) ? η_m/ν : ηj_p*ν)
    return real.(B + B_spin * sdot)
end

function site_field_eta(mc, pos)
    B = zeros(3)
    for disp in disps
        posj = pos .+ disp
        B += bond_field_eta(mc, Dimer(pos, posj))
    end
    return B
end

## Spin sweep functions ##

function biasfield(_::WignerMC{AlgType, Nothing}, _, _) where AlgType
    return [0,0,0]
end
biasfield(mc::WignerMC, pos, B) = B * mc.bias(pos...)

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

function site_energy_s(mc::WignerMC, pos, s, B)
    E = 0.0
    for disp in disps
        posj = pos + disp
        sdot = s ⋅ mc.spins[posj...] / 4
        E += bond_energy_s(mc, Dimer(pos, posj), sdot)
    end
    return E - biasfield(mc, pos, B) ⋅ s
end

"""
    bond_field_s(mc::WignerMC, pos, posj)

Calculate spin field on `pos` due to `posj`
"""
bond_field_s(mc::WignerMC, pos, posj) = ssfactor(mc, Dimer(pos, posj)) * mc.spins[posj...] / 2

function site_field_s(mc::WignerMC, pos)
    B = [0,0,0]
    for disp in disps
        posj = pos + disp
        B += bond_field_s(mc, pos, posj)
    end
    return B
end
function site_field_s(mc::WignerMC, pos, B)
    return site_field_s(mc, pos) - biasfield(mc, pos, B)
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
    return E
end
half_energy(mc, pos) = half_energy_nobias(mc, pos)
function half_energy(mc::WignerMC, pos)
    half_energy_nobias(mc, pos) - biasfield(mc, pos, mc.B) ⋅ mc.spins[pos...] / 2
end

function total_energy(mc)
    E = 0.0
    for I in eachindex(IndexCartesian(), mc.spins)
        E += half_energy(mc, Tuple(I))
    end
    return E
end
