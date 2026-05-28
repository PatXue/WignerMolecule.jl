function bond_energy(mc::EtaMC, s::SpinVector, sj::SpinVector, ν)
    # Coupling energies
    Jzz = mc.params.Jzz
    Jp = mc.params.Jp

    s /= 2
    sj /= 2

    # η raising and lowering operators
    s_p = s[1] + 1.0im*s[2]
    sj_m = sj[1] - 1.0im*sj[2]
    sj_p = sj[1] + 1.0im*sj[2]

    E = 0.0
    E += Jzz * s[3] * sj[3]
    E += real(J * s_p * sj_m)
    E += 2Jp * real(ν * s_p * sj_p)

    return -E
end

