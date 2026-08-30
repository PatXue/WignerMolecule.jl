ismonomer(pos, mc::DimerMC) = mod_equiv(mc.spins[pos...], pos, mc)
indimer(pos, d::Dimer, mc::DimerMC) = mod_equiv(pos, d.pos, mc) || mod_equiv(pos, d.posj, mc)

# Get dimer pos is part of (assumes pos is not a monomer)
getdimer(pos, mc::DimerMC) = Dimer(pos, mc.spins[pos...])
"""
    arestacked(d::Dimer, dj::Dimer, mc)

Check if double dimer move is possible (two parallel stacked dimers). Dimers must be oriented
"""
arestacked(d::Dimer, dj::Dimer, mc) = areneighbors(d.pos, dj.pos, mc) && areneighbors(d.posj, dj.posj, mc)

function add_dimer!(d::Dimer, mc::DimerMC)
    mc.spins[d.pos...] = d.posj
    mc.spins[d.posj...] = d.pos
end

# Monomer BitSet handling functions
function postoint(pos, mc::DimerMC)
    Lx = size(mc.spins, 1)
    pos = mod.(pos, size(mc.spins))
    return pos[1] + Lx * pos[2]
end
function inttopos(n, mc::DimerMC)
    Lx = size(mc.spins, 1)
    return SVector{2,Int}(n % Lx, fld(n, Lx))
end
function addmonomer!(pos, s, mc::DimerMC)
    mc.monospins[pos...] = s
    push!(mc.monomers, postoint(pos, mc))
end
function delmonomer!(pos, mc::DimerMC)
    mc.monospins[pos...] = SVector(0,0,0)
    pop!(mc.monomers, postoint(pos, mc))
end

randmonomer(mc::DimerMC, rng=default_rng()) = inttopos(rand(rng, mc.monomers), mc)

# Perform Fourier transform on MC, updating preallocated ηks matrix
function update_fourier!(mc::DimerMC)
    Lx, Ly = size(mc.spins)
    for (x, y) in Iterators.product(1:Lx, 1:Ly)
        disp = mc.spins[x,y] - [x,y]
        mc.sks[x,y,1:3] .= [mod_equiv(disp, oriented_disps[i], mc) for i in 1:3]
        mc.sks[x,y,4] = mod_equiv(disp, [0,0], mc)
    end
    fft!(mc.sks.array, (1, 2))
    mc.sks ./= length(mc.spins)

    for i in 1:3
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.ηks.array, (1, 2))
    mc.ηks ./= length(mc.ηs)
    return nothing
end

function validate_mc(mc::DimerMC)
    Lx, _ = size(mc.spins)
    for I in eachindex(mc.spins)
        x,y = Tuple(I)
        pos = SVector(x, y)
        posj = mc.spins[pos...]
        @assert mod_equiv(pos, mc.spins[posj...], mc)
        if ismonomer(pos, mc)
            @assert mod_equiv(pos, posj, mc)
            @assert (x + Lx * y) in mc.monomers
            @assert norm(mc.monospins[pos...]) ≈ 1.0
        end
    end
    return true
end
