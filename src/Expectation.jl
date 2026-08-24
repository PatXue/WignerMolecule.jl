module Expectations

export Expectation, combine, addsample

struct Expectation
    val::Float64    # Value of measured expectation
    n::Int          # Number of samples expectations calculated from
end

"""
    combine(e1::Expectation, e2::Expectation)

Combine two measured expectations (carrying out error propagation assuming
independence)
"""
function combine(e1::Expectation, e2::Expectation)
    n1 = e1.n
    n2 = e2.n
    val = (n1 * e1.val + n2 * e2.val) / (n1 + n2)
    return Expectation(val, n1+n2)
end

function addsample(e::Expectation, x)
    if e.n == 0
        return Expectation(x, 1)
    end
    val = e.val
    n = e.n
    new_val = (n*val + x) / (n+1)
    return Expectation(new_val, n+1)
end

end