function bond_energy(mc::DimerMC, I, J, ν)
    # Couplings
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

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
    E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
    E_spin += 2*J_EMEP_SS *     η_m * ηj_p
    E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
    E_spin += J_SS
    E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)

    return real(E_spin * (s⋅sj) + E_η)
end

