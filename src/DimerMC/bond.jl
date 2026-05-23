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
function getpartner(mc::DimerMC, (x, y))
    b = mc.spins[x, y]
    return (x, y) .+ bondtodisp[b]
end
