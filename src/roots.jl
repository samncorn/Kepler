function bisect_step(xl, xh)
    return 0.5(xl + xh)
end

function lmm12_step(x1, x2, f1, f2, df1, df2)
    # construct polynomial approximation of the inverse
    # cubic in this case
    p = SVector{4}((x1, x2, 1/df1, 1/df2))
    M = transpose(SMatrix{4, 4}(
        1.0,  f1, f1^2,   f1^3,
        1.0,  f2, f2^2,   f2^3,
        0.0, 1.0, 2*f1, 3*f1^2,
        0.0, 1.0, 2*f2, 3*f2^2,
    ))
    # interpolate the inverse
    # f^(-1)(x) = b0 + b_i x^i, so
    # f^(-1)(0) = b0
    b0, _, _, _ = pinv(M)*p

    if min(x1, x2) < b0 < max(x1, x2) # need additional check to safeguard against slow convergence
        return b0
    else
        return 0.5(x1 + x2)
    end
end

function check_step(x, xl, xh; tol = 1e-14)
    return x == xl || x == xh || abs(xl - xh) < tol
end