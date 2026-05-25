# Helper functions for handling Bonds

struct Dimer
    x::Int,
    y::Int,
    b::Bond
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
getpartner(d::Dimer) = (d.x, d.y) .+ bondtodisp(d.b)

# Perform a rotation of x hat to lie along given bond r
rotate(b::Bond, r::Bond) = Bond((Int(b) + Int(r)) % 6)
function rotate((x, y), r::Bond)
    x .* bondtodisp[r] .+ y .* bondtodisp[rotate(a2, r)]
end
rotate(d::Dimer, r::Bond) = Dimer(rotate((d.x,d.y), r)..., rotate(d.b, r))
