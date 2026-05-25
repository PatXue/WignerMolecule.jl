# Helper functions for handling Bonds

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

# Rotate x hat to the given bond
rotate(b::Bond, r::Bond) = Bond((Int(b) + Int(r)) % 6)
function rotate((x, y), r::Bond)
    x .* bondtodisp[r] .+ y .* bondtodisp[rotate(a2, r)]
end
