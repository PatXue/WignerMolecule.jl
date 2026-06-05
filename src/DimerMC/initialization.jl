# Initialize the DimerMC spins/etas
const unitcell = [(0,1) (0,-1); (1,-1) (-1,1)]
function init_vbs_s!(mc::DimerMC)
    for I in eachindex(IndexCartesian(), mc.spins)
        x, y = Tuple(I)
        # Shift x,y to be inside unit cell
        while y > 2
            x += 1
            y -= 2
        end
        x = mod1(x, 2)
        mc.spins[I] = SVector((Tuple(I) .+ unitcell[x,y])...)
    end
end

function Carlo.init!(mc::DimerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_vbs_s!(mc)
    rand!(mc.ηs)
    update_fourier!(mc)
    return nothing
end
