# Helper functions for handling Bonds

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

# Shift position by v
shift(d::Dimer, v) = Dimer(d.pos .+ v, d.posj .+ v)

# Reflect position across x-axis
const reflect1_mat = SMatrix{2,2}(1, 0, 1, -1)
reflect1(d::Dimer) = Dimer(reflect1_mat * d.pos, reflect1_mat * d.posj)
# Reflect position across line 30 deg above x-axis
const reflect2_mat = SMatrix{2,2}(0, 1, 1, 0)
reflect2(d::Dimer) = Dimer(reflect2_mat * d.pos, reflect2_mat * d.posj)

# Rotation matrices in increments of 60 degrees
const rotmats = Tuple(SMatrix{2,2}(0, 1, -1, 1)^i for i in 0:5)

# Rotate around (0,0) by r*60 degrees
rotate(v, r) = rotmats[r+1] * v
rotate(d::Dimer, r) = Dimer(rotate(d.pos, r), rotate(d.posj, r))
invrotate(x, r) = rotate(x, mod(6-r, 6))

# Return dimers that d conflicts with in mc
function collisions(mc::DimerMC, d::Dimer)
    res = []
    if mod_equiv(mc.spins[d.pos...], d.posj, mc)
        return res
    end
    if !mc.visited[d.pos...]
        push!(res, Dimer(d.pos, mc.spins[d.pos...]))
    end
    if !mc.visited[d.posj...]
        push!(res, Dimer(d.posj, mc.spins[d.posj...]))
    end
    return res
end

# Check and flip dimer to lie along a1, a2, or a3
const oriented_disps = (SVector(1,0), SVector(-1,1), SVector(0,-1))
function orientdimer(d::Dimer, mc::DimerMC)
    Lx, Ly == size(mc.spins)
    disp = mod1(d.posj - d.pos, mc)
    if any([mod_equiv(disp, a, mc) for a in oriented_disps])
        return d
    else
        return Dimer(d.posj, d.pos)
    end
end
# Get the ν coupling factor for a dimer (assuming dimer oriented)
function getν(d::Dimer, mc::DimerMC)
    Lx, Ly == size(mc.spins)
    disp = mod1(d.posj - d.pos, mc)
    if mod_equiv(disp, (1,0), mc)
        return 1
    elseif mod_equiv(disp, (-1,1), mc)
        return ω
    elseif mod_equiv(disp, (0,-1), mc)
        return ω^2
    else
        throw(ArgumentError("Dimer $d not oriented correctly"))
    end
end
