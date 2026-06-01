# Helper functions for handling Bonds

struct Dimer
    pos::SVector{2, Int}
    posj::SVector{2, Int}
end

# Apply periodic boundaries to positions
Base.mod1(pos::SVector{2,Int}, mc::DimerMC) = mod1.(pos, size(mc.spins))
function Base.mod1(d::Dimer, mc::DimerMC)
    dims = size(mc.spins)
    return Dimer(mod1.(d.pos, dims), mod1.(d.posj, dims))
end

# Shift position by v
shift(d::Dimer, v) = Dimer(d.pos .+ v, d.posj .+ v)

# Reflect position across x-axis
const reflect1_mat = SMatrix{2,2}(1, 0, 1, -1)
reflect1(d::Dimer) = Dimer(reflect1_mat * d.pos, reflect1_mat * d.posj)
# Reflect position across line 30 deg above x-axis
const reflect2_mat = SMatrix{2,2}(0, 1, 1, 0)
reflect2(d::Dimer) = Dimer(reflect2_mat * d.pos, reflect2_mat * d.posj)

# Rotation matrices in increments of 60 degrees
const rotmats = @SVector [SMatrix{2,2}(0, 1, -1, 1)^i for i in 0:5]

# Rotate around (0,0) by r*60 degrees
rotate(v, r) = rotmats[r+1] * v
rotate(d::Dimer, r::Bond) = Dimer(rotate(d.pos, r), rotate(d.posj, r))
invrotate(x, r::Bond) = rotate(x, reflect1(r))

# Return dimers that d conflicts with in mc
function collisions(mc::DimerMC, d::Dimer)
    res = []
    if mc.spins[d.pos...] == mod1(d.posj, mc)
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

# Check and flip dimer to lie along a1, a3, or a5
function orientdimer(d::Dimer)
    if disptobond[d.posj - d.pos] ∈ (a1, a3, a5)
        return d
    else
        return Dimer(d.posj, d.pos)
    end
end
# Get the ν coupling factor for a dimer
function getν(d::Dimer)
    bond = disptobond[d.posj - d.pos]
    if bond ∈ (a1, a4)
        return 1
    elseif bond ∈ (a2, a5)
        return ω
    elseif bond ∈ (a3, a6)
        return ω^2
    end
end
