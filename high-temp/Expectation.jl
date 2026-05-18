module Expectations

export Expectation, combine, addsample

struct Expectation
    val::Float64    # Value of measured expectation
    err::Float64    # Std. dev in measured expectation
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
    err = sqrt((n1 * e1.err)^2 + (n2 * e2.err)^2) / (n1 + n2)
    return Expectation(val, err, n1+n2)
end

function addsample(e::Expectation, x)
    if e.n == 0
        return Expectation(x, 0.0, 1)
    end
    val = e.val
    err = e.err
    n = e.n
    new_val = (n*val + x) / (n+1)
    old_sqs = (n-1) * err^2 + val^2  # Avg of squares for old samples
    new_sqs = (n*old_sqs + x^2) / (n+1)
    return Expectation(new_val, sqrt((new_sqs - new_val^2) / n), n+1)
end

end