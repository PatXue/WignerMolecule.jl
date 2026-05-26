# Helper functions for handling Bonds

struct Dimer
    pos::Tuple{Int, Int}
    posj::Tuple{Int, Int}
end

# Bond type to index displacement
const bondtodisp = Dict(
    a1 => (1,0),
    a2 => (0,1),
    a3 => (-1,1),
    a4 => (-1,0),
    a5 => (0,-1),
    a6 => (1,-1)
)

# Get position index of (x,y)'s entanglement partner
getpartner(mc::DimerMC, (x, y)) = (x, y) .+ bondtodisp[mc.spins[x, y]]

# Perform a rotation of x hat to lie along given bond r
rotate(b::Bond, r::Bond) = Bond((Int(b) + Int(r)) % 6)
function rotate((x, y), r::Bond)
    x .* bondtodisp[r] .+ y .* bondtodisp[rotate(a2, r)]
end
rotate(d::Dimer, r::Bond) = Dimer(rotate(d.pos, r), rotate(d.posj, r))

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
