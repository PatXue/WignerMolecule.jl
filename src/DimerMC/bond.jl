# Helper functions for handling indexing and dimers
struct Dimer
    pos::SVector{2, Int}
    posj::SVector{2, Int}
end

# Apply periodic boundaries to positions
Base.mod1(pos, mc::DimerMC) = mod1.(pos, size(mc.spins))
function Base.mod1(d::Dimer, mc::DimerMC)
    dims = size(mc.spins)
    return Dimer(mod1.(d.pos, dims), mod1.(d.posj, dims))
end
mod_equiv(pos, posj, mc::DimerMC) = all((pos .- posj) .% size(mc.spins) .== (0,0))

const disps = (SVector(1,0), SVector(-1,1), SVector(0,-1), SVector(-1,0), SVector(1,-1), SVector(0,1))
const oriented_disps = (SVector(1,0), SVector(-1,1), SVector(0,-1))
# Check and flip dimer to lie along a1, a2, or a3
function orientdimer(d::Dimer, mc::DimerMC)
    disp = d.posj - d.pos
    if any([mod_equiv(disp, a, mc) for a in oriented_disps])
        return d
    else
        return Dimer(d.posj, d.pos)
    end
end
# Get the ν coupling factor for a dimer (assuming dimer oriented)
function getν(d::Dimer, mc::DimerMC)
    disp = d.posj - d.pos
    getν(disp, mc)
end
function getν(disp, mc::DimerMC)
    if mod_equiv(disp, (1,0), mc)
        return 1
    elseif mod_equiv(disp, (-1,1), mc)
        return ω
    elseif mod_equiv(disp, (0,-1), mc)
        return ω^2
    else
        throw(ArgumentError("Displacement $disp invalid"))
    end
end

ismonomer(pos, mc::DimerMC) = mod_equiv(mc.spins[pos...], pos, mc)
indimer(pos, d::Dimer, mc::DimerMC) = mod_equiv(pos, d.pos, mc) || mod_equiv(pos, d.posj, mc)

# Monomer BitSet handling functions
function addmonomer!(pos, s::SpinVector, mc::DimerMC)
    Lx = size(mc.spins, 1)
    pos = mod1.(pos, size(mc.spins))
    mc.spins[pos...] = pos
    mc.monospins[pos...] = s
    push!(mc.monomers, pos[1] + Lx * pos[2])
end
function delmonomer!(pos, mc::DimerMC)
    Lx = size(mc.spins, 1)
    pos = mod1.(pos, size(mc.spins))
    pop!(mc.monomers, pos[1] + Lx * pos[2])
end

function randmonomer(mc::DimerMC, rng=default_rng())
    Lx = size(mc.spins, 1)
    n = rand(rng, mc.monomers)
    return (div(n, Lx), n % Lx)
end
