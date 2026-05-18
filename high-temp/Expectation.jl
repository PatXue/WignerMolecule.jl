module Expectations

export Expectation, combine

struct Expectation
    val::Float64    # Value of measured expectation
    err::Float64    # Std. dev in measured expectation
    n::UInt64       # Number of samples expectations calculated from
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
    err = sqrt((n1 * e1.err)^2 + (n2 * e2.err)^2) / (n1 + n2)
    return Expectation(val, err, n1+n2)
end

end