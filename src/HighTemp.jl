# Sampling functions for high-T expectations
module HighTemp

export randomize!

# Randomize spins in MC
function randomize!(mc::WignerMC)
    rand!(mc.spins)
    rand!(mc.ηs)
end

# Sample values of functions in fs over n runs on random spins
function sample!(mc::WignerMC, n, fs::Vararg{Function, N}) where {N}
    means = fill(0.0, N)
    errs = fill(0.0, N)
    for _ in 1:n
        randomize!(mc)
        vals = [f(mc) for f in fs]
        @. means = means + vals
        @. errs = errs + vals^2
    end
    @. return (means / n, sqrt((errs - means) / (n-1)))
end

end
