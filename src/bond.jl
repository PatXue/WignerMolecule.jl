# Helper functions for handling indexing and bonds
struct Dimer
    pos::SVector{2, Int}
    posj::SVector{2, Int}
end

# Apply periodic boundaries to positions
Base.mod1(pos, mc) = mod1.(pos, size(mc.spins))
function Base.mod1(d::Dimer, mc)
    dims = size(mc.spins)
    return Dimer(mod1.(d.pos, dims), mod1.(d.posj, dims))
end
mod_equiv(pos, posj, mc) = all((pos .- posj) .% size(mc.spins) .== (0,0))

const disps = (SVector(1,0), SVector(-1,1), SVector(0,-1), SVector(-1,0), SVector(1,-1), SVector(0,1))
const oriented_disps = (SVector(1,0), SVector(-1,1), SVector(0,-1))

isoriented(d::Dimer, mc) = any([mod_equiv(d.posj - d.pos, a, mc) for a in oriented_disps])
# Check and flip dimer to lie along a1, a2, or a3
function orientdimer(d::Dimer, mc)
    if isoriented(d, mc)
        return d
    else
        return Dimer(d.posj, d.pos)
    end
end

# Get the ν coupling factor for a dimer (assuming dimer oriented)
function getν(d::Dimer, mc)
    disp = d.posj - d.pos
    getν(disp, mc)
end
function getν(disp, mc)
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
