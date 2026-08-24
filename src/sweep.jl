function metropolisacc(ΔE, T; rng, fug=1.0)
    prob = fug * exp(-ΔE/T)
    return prob >= 1.0 || rand(rng) <= prob
end

"""
    heatbath_spin(Bvec, T; rng)

Given the local field Bvec on a site, sample a spin from the local Boltzmann distribution
"""
function heatbath_spin(Bvec, T; rng)
    B = norm(Bvec)
    if isapprox(B/T, 0.0, atol=1e-5)
        return rand(rng, SpinVector)
    end

    # Generate parallel component
    Bhat = Bvec / B
    c = B / 2T
    cosθ = -log(exp(c) - 2*rand(rng)*sinh(c)) / c
    s = Bhat * cosθ
    # Find a basis of B's orthogonal subspace
    u = Bhat × SVector(1,0,0)
    if isapprox(norm(u), 0.0, atol=1e-4)
        u = Bhat × SVector(0,1,0)
    end
    u /= norm(u)
    v = Bhat × u
    # Generate orthogonal component
    ϕ = 2pi * rand(rng)
    s += sqrt(1-norm(s)^2) * (u*cos(ϕ) + v*sin(ϕ))

    @assert norm(s) ≈ 1.0
    return s
end

function flip_eta!(mc::WignerMC, T, rng)
    Lx, Ly = size(mc.ηs)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

    old_E = site_energy_eta(mc, pos, mc.ηs[pos...])
    new_η = rand(rng, SpinVector)
    new_E = site_energy_eta(mc, pos, new_η)
    ΔE = new_E - old_E

    if metropolisacc(ΔE, T; rng)
        mc.ηs[pos...] = new_η
    end
end

function flip_eta!(mc::WignerMC{:Heatbath}, T, rng)
    Lx, Ly = size(mc.ηs)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    mc.ηs[pos...] = heatbath_spin(site_field_eta(mc, pos), T; rng)
end

function flip_spin!(mc::WignerMC, T, B, rng)
    Lx, Ly = size(mc.ηs)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))

    old_E = site_energy_s(mc, pos, mc.spins[pos...], B)
    new_s = rand(rng, SpinVector)
    new_E = site_energy_s(mc, pos, new_s, B)
    ΔE = new_E - old_E

    if metropolisacc(ΔE, T; rng)
        mc.spins[pos...] = new_s
    end
end

function flip_spin!(mc::WignerMC{:Heatbath}, T, B, rng)
    Lx, Ly = size(mc.spins)
    pos = SVector(rand(rng, 1:Lx), rand(rng, 1:Ly))
    mc.spins[pos...] = heatbath_spin(site_field_s(mc, pos, B), T; rng)
end

function Carlo.sweep!(mc::WignerMC, ctx::Carlo.MCContext)
    rng = ctx.rng
    T = calc_temp(mc, ctx)
    for _ in 1:length(mc.spins)
        flip_eta!(mc, T, rng)
        flip_spin!(mc, T, calc_B(mc, ctx), rng)
    end
    return nothing
end


function cluster_spin!(mc, T, rng)
    Lx, Ly = size(mc.spins)
    nhat = rand(SpinVector, rng)
    pos = SVector{2,Int}(rand(rng, 1:Lx), rand(rng, 1:Ly))

    cluster = Set{SVector{2,Int}}()
    newposns = Set{SVector{2,Int}}()
    push!(cluster, pos)
    push!(newposns, pos)

    while length(newposns) > 0
        pos = pop!(newposns)
        s = mc.spins[pos...] / 2
        for a in disps
            posj = mod1.(pos + a, (Lx, Ly))
            if posj in cluster
                continue
            end
            sj = mc.spins[posj...] / 2
            prob = 1 - exp(min(0, 2 * ssfactor(mc, Dimer(pos, posj)) * (nhat ⋅ s) * (nhat ⋅ sj) / T))
            if rand(rng) < prob
                push!(cluster, posj)
                push!(newposns, posj)
            end
        end
    end

    for pos in cluster
        mc.spins[pos...] -= 2 * (nhat ⋅ mc.spins[pos...]) * nhat
    end
    return length(cluster)
end

function Carlo.sweep!(mc::WignerMC{:Cluster}, ctx::Carlo.MCContext)
    rng = ctx.rng
    T = calc_temp(mc, ctx)
    for _ in 1:length(mc.spins)
        flip_eta!(mc, T, rng)
        cluster_spin!(mc, T, rng)
    end
    return nothing
end
