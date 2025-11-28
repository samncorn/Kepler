function bisect_step(xl, xh)
    return 0.5(xl + xh)
end

function lmm12_step(x1, x2, f1, f2, df1, df2)
    # construct polynomial approximation of the inverse
    # cubic in this case
    # p = SVector{4}((x1, x2, 1/df1, 1/df2))
    # M = transpose(SMatrix{4, 4}(
    #     1.0,  f1, f1^2,   f1^3,
    #     1.0,  f2, f2^2,   f2^3,
    #     0.0, 1.0, 2*f1, 3*f1^2,
    #     0.0, 1.0, 2*f2, 3*f2^2,
    # ))
    # # interpolate the inverse
    # # f^(-1)(x) = b0 + b_i x^i, so
    # # f^(-1)(0) = b0
    # b0, _, _, _ = pinv(M)*p
    # return b0

    t   = -f1/(f2 - f1)
    h00 = (1 + 2t)*(1 - t)^2
    h10 = t*(1-t)^2
    h01 = (3 - 2t)*t^2
    h11 = (t - 1)*t^2

    return h00*x1 + h10*(f2-f1)/df1 + h01*x2 + h11*(f2-f1)/df2
end

# function cubic_hermite_spline(x, x1, x2, y1, y2)


# function brent_lmm12_step(xa, del1, del2, x1, x2, f1, f2, df1, df2, previous_bisect; tol = 1e-15)
#     if previous_bisect && abs(del1) < tol
#         return 0.5(xa + x1), true
#     end

#     if !previous_bisect && abs(del2) < tol
#         return 0.5(xa + x1), true
#     end

#     s = lmm12_step(x1, x2, f1, f2, df1, df2)
#     if previous_bisect && 
#         return 0.5(xa + x1), true
#     end

#     if !previous_bisect && !((3xa + x1)/4 < s < x1)
# end


function check_step(x, xl, xh; tol = 1e-14)
    return x == xl || x == xh || abs(xl - xh) < tol
end