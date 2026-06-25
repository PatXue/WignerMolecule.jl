function bond_energy(mc::EtaMC, s::SpinVector, sj::SpinVector, ν)
    # Coupling energies
    Jzz = mc.params.Jzz
    Jp = mc.params.Jp

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

# Calculate the energy at a lattice site (x, y) if it had spin s
function energy(mc::EtaMC, s::SpinVector, x, y)
    E = 0.0
    # Nearest neighbor lattice positions along a1,a2,a3
    nns = ((x+1, y), (x-1, y+1), (x, y-1))
    for j in eachindex(nns)
        ν = ω^(j-1)
        nn = nns[j]
        sj = mc.spins[nn...]
        E += bond_energy(mc, s, sj, ν)
    end

    # Nearest neighbor lattice positions along -a1,-a2,-a3
    nns = ((x-1, y), (x+1, y-1), (x, y+1))
    for j in eachindex(nns)
        ν = ω^(j-1)
        nn = nns[j]
        sj = mc.spins[nn...]
        E += bond_energy(mc, sj, s, ν)
    end
    return E - mc.B * s[3]
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::EtaMC, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x-1, y+1), (x, y-1))
    E = 0.0
    s = mc.spins[x, y]
    for j in eachindex(nns)
        ν = ω^(j-1)
        nn = nns[j]
        sj = mc.spins[nn...]
        E += bond_energy(mc, s, sj, ν)
    end
    return E - mc.B * s[3]
end

# Calculate the total energy of MC
function total_energy(mc::EtaMC)
    tot_energy = 0.0
    for I in eachindex(mc.spins)
        x, y = Tuple(I)
        tot_energy += half_energy(mc, x, y)
    end
    return tot_energy
end
