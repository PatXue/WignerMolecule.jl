# Helper functions for calculating system energy

function bond_energy(mc::DimerMC, paired, η, ηj, ν)
    # Couplings
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    η /= 2
    ηj /= 2

    # η raising and lowering operators
    η_m = η[1] + 1.0im*η[2]
    ηj_p = ηj[1] - 1.0im*ηj[2]
    ηj_m = ηj[1] + 1.0im*ηj[2]

    E_spin = 0.0 + 0.0im
    E_η = 0.0 + 0.0im

    # η-only energy
    E_η +=   J_EzEz *     η[3] * ηj[3]
    E_η += 2*J_EMEP *     η_m * ηj_p
    E_η += 2*J_EMEM * ν * η_m * ηj_m

    # η-S energy
    if paired
        E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
        E_spin += 2*J_EMEP_SS *     η_m * ηj_p
        E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
        E_spin += J_SS
        E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)
        E_spin *= -3/4
    end

    return real(E_spin + E_η)
end

function site_energy(mc::DimerMC, η, pos)
    E = 0.0
    for j in 1:3
        ν = ω^(j-1)
        disp = oriented_disps[j]

        posj = pos .+ disp
        ηj = mc.ηs[posj...]
        paired = mod_equiv(mc.spins[pos...], posj, mc)
        E += bond_energy(mc, paired, η, ηj, ν)

        posj = pos .- disp
        ηj = mc.ηs[posj...]
        paired = mod_equiv(mc.spins[pos...], posj, mc)
        E += bond_energy(mc, paired, ηj, η, ν)
    end
    return E
end


# Energy from spin-orbit coupling if d were entangled
function bond_energy(mc::DimerMC, d::Dimer)
    # Couplings
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS

    d = orientdimer(d, mc)
    η = mc.ηs[d.pos...] / 2
    ηj = mc.ηs[d.posj...] / 2
    ν = getν(d, mc)

    # η raising and lowering operators
    η_m = η[1] + 1.0im*η[2]
    ηj_p = ηj[1] - 1.0im*ηj[2]
    ηj_m = ηj[1] + 1.0im*ηj[2]

    E = 0.0 + 0.0im
    E +=   J_EzEz_SS *     η[3] * ηj[3]
    E += 2*J_EMEP_SS *     η_m * ηj_p
    E += 2*J_EMEM_SS * ν * η_m * ηj_m
    E += J_SS
    E += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)
    E *= -3/4

    return real(E)
end
