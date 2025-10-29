const ω::ComplexF64 = exp(im * 2π/3)
# Calculate the energy at a lattice site (x, y) if it had spin s and
# pseudospin η
function energy(mc::WignerMC, s::SpinVector, η::SpinVector, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1), (x-1, y), (x-1, y+1), (x, y+1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

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

        E += real(E_spin * (s⋅sj) + E_η)
    end

    return E
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::WignerMC, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    s = mc.spins[x, y]
    η = mc.ηs[x, y]
    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

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

        E += real(E_spin * (s⋅sj) + E_η)
    end

    return E
end

# Calculate the energy difference at a lattice site (x, y) if the spin changed
# by s_diff
function s_energydiff(mc::WignerMC, s_diff::SpinVector, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1), (x-1, y), (x-1, y+1), (x, y+1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS

    # neighbor spins weighted by η couplings
    weighted_sj = [0, 0, 0]
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

        # η raising and lowering operators
        η_m = η[1] + 1.0im*η[2]
        ηj_p = ηj[1] - 1.0im*ηj[2]
        ηj_m = ηj[1] + 1.0im*ηj[2]

        # η-S energy
        E_spin = 0.0 + 0.0im
        E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
        E_spin += 2*J_EMEP_SS *     η_m * ηj_p
        E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
        E_spin += J_SS
        E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)

        weighted_sj .+= E_spin .* sj
    end

    return s_diff ⋅ weighted_sj
end

# Calculate the energy difference at a lattice site (x, y) if the η changed
# by η_diff
function η_energydiff(mc::WignerMC, η_diff::SpinVector, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1), (x-1, y), (x-1, y+1), (x, y+1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

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

        E += real(E_spin * (s⋅sj) + E_η)
    end

    return E
end
