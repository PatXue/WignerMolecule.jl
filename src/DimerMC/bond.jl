# Helper functions for handling Bonds

struct Dimer
    pos::SVector{2, Int}
    posj::SVector{2, Int}
end

# Bond type to index displacement
const bondtodisp = Dict(
    a1 => SVector(1,0),
    a2 => SVector(0,1),
    a3 => SVector(-1,1),
    a4 => SVector(-1,0),
    a5 => SVector(0,-1),
    a6 => SVector(1,-1)
)

# Get position index of v = (x,y)'s entanglement partner in mc
getpartner(mc::DimerMC, v) = v .+ bondtodisp[mc.spins[v[1], v[2]]]

# Bond type to rotation matrix
const bondtorot = Dict(
    a1 => SMatrix{2,2}(1, 0, 0, 1),
    a2 => SMatrix{2,2}(0, 1, -1, 1),
    a3 => SMatrix{2,2}(-1, 1, -1, 0),
    a4 => SMatrix{2,2}(-1, 0, 0, -1),
    a5 => SMatrix{2,2}(0, -1, 1, -1),
    a6 => SMatrix{2,2}(1, -1, 1, 0)
)

# Perform a rotation of x hat to lie along given bond r
rotate(b::Bond, r::Bond) = Bond((Int(b) + Int(r)) % 6)
rotate(v, r::Bond) = bondtorot[r] * v
rotate(d::Dimer, r::Bond) = Dimer(rotate(d.pos, r), rotate(d.posj, r))

# Reflect position across x-axis
reflect1((x, y),) = x .* (1,0) .+ y .* (1,-1)
# Reflect position across line 30 deg above x-axis
reflect2((x, y),) = x .* (0,1) .+ y .* (1,0)

# Return dimers that d conflicts with in mc
function collisions(mc::DimerMC, d::Dimer)
    res = []
    if getpartner(mc, d.pos) == d.posj
        return res
    end
    if mc.spins[d.pos...] != none
        push!(res, Dimer(d.pos, getpartner(mc, d.pos)))
    end
    if mc.spins[d.posj...] != none
        push!(res, Dimer(d.posj, getpartner(mc, d.posj)))
    end
    return res
end
