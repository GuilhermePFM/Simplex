"""
Return a copy of `lp` with all right-hand side values non-negative.
"""
function make_b_nonnegative(lp::LinearProgram) :: LinearProgram
    mask = lp.b .< 0
    any(mask) || return lp
    A2 = copy(lp.A)
    b2 = copy(lp.b)
    A2[mask, :] .*= -1
    b2[mask]    .*= -1
    LinearProgram(A2, b2, lp.c)
end
