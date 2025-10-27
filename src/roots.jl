
"""
f needs to return the value and first derivative
f also needs to be reasonably well behaved on the bracket interval (possibly strictly increasing or decreasing)

This is guarunteed to converge if a root is bracketed, with the usual rootfinding caveats (may be slow in some cases, 
might find multiple roots, etc). 
"""
function brent_newton(f, br; max_iters = 1000, tol = 1e-15, btol = 1e-8)
    # choose the contrapoint (b)
    a, b   = br
    fa, va = f(a)
    fb, vb = f(b)
    @assert sign(fa) != sign(fb) "root not bracketed"

    if abs(fa) > abs(fb)
        # swap them
        a, b   = (b, a)
        fa, fb = (fb, fa)
        va, vb = (vb, va)
    end

    if fa == 0.0
        return a
    end

    i   = 0
    # xp  = Inf # x n-1
    # xpp = Inf # x n-2
    dx1 = abs(a - b)
    dx2 = dx1
    # max_iters = 2(ceil(Int, log2(dx1/tol))^2)
    max_iters = Inf
    # @debug "estimated maximum $(max_iters) iterations to converge"
    # bs  = false # whether to use bisection
    bs = false # whatever the last step used
    while i < max_iters && abs(fa) > tol
        i += 1
        # x = bs ? 0.5(x1 + x2) : interpolate(0.0, Cubic_Hermite_Spline(y1, y2, x1, x2, 1/v1, 1/v2))
        bs = false
        x  = interpolate(0.0, Cubic_Hermite_Spline(fa, fb, a, b, 1/va, 1/vb))
        dx = abs(x - a)

        # check if we are making adquate progress
        if !((bs && btol < dx1 && 2dx < dx1) || (!bs && btol < dx2 && 2dx < dx2))
            # bisect
            x  = 0.5(a + b)
            dx = abs(x - a)
            bs = true
        end

        # compute the new value
        fx, vx = f(x)

        # set up the new interval
        dx2 = dx1
        dx1 = dx

        # find which side the root is on, choose new contrapoint
        if sign(fx) == sign(fa)
            if abs(fx) < abs(fb)
                a  = x
                fa = fx
                va = vx
            else
                a  = b
                fa = fb
                va = vb
                b  = x
                fb = fx
                va = vx
            end
        elseif sign(fx) == sign(fb)
            if abs(fx) < abs(fa)
                b  = a
                fb = fa
                vb = va
                a  = x
                fa = fx
                va = vx
            else
                b  = x
                fb = fx
                va = vx
            end
        else
            # we have found the root (sign (0) == 0)
            return x
        end
    end
    i >= max_iters && throw("exceeded max iterations $(max_iters)")
    return a
end

function brent_newton_step(br_x, br_f, br_df, dx1, dx2, bs; btol = 1e-8)
    @assert sign(br_f[1]) != sign(br_f[2]) 
    x = 0.5(sum(br_x))
    if br_df[1] != 0 && br_df[2] != 0
        # interpolate x(f)
        x  = interpolate(0.0, Cubic_Hermite_Spline(br_f..., br_x..., 1/br_df[1], 1/br_df[2]))
        dx = br_f[1] < br_f[2] ? abs(x - br_x[1]) : abs(x - br_x[2])
        x, dx
    else
        # bisect
        x   = 0.5sum(br_x)
        return x, true
    end

    # check progress
    if !((bs && btol < dx1 && 2dx < dx1) || (!bs && btol < dx2 && 2dx < dx2))
        # bisect
        x   = 0.5sum(br_x)
        return x, true
    end

    return x, false
end

struct Cubic_Hermite_Spline{T}
    endpoints   ::Tuple{T, T}
    coefficients::Tuple{T, T, T, T}
end

function Cubic_Hermite_Spline(x1, x2, y1, y2, v1, v2)
    if x2 < x1
        return Cubic_Hermite_Spline(x2, x1, y2, y1, v2, v1)
    end
    v1 = v1 * (x2 - x1)
    v2 = v2 * (x2 - x1)
    coeff = (
        y1,
        v1,
        3(y2 - y1) - 2v1 - v2,
        -2(y2 - y1) + v1 + v2
    )
    return Cubic_Hermite_Spline((x1, x2), coeff)
end

function interpolate(x, spline)
    x1, x2         = spline.endpoints
    a1, a2, a3, a4 = spline.coefficients
    t = (x - x1)/(x2 - x1)
    return a1 + t*(a2 + t*(a3 + t*a4))
end

""" given computed values for f and it's derivative (df) at 2 points, computes the next step by inverse 
interpolation. Analogous to a Newton-Raphson step, but with two history points.
"""
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
    return b0
end

"""given computed values for f and it's derivative (df) at 3 points, computes the next step by inverse 
interpolation. Analogous to a Newton-Raphson step, but with three history points.
"""
function lmm13_step(x1, x2, x3, f1, f2, f3, df1, df2, df3)
    # construct (quintic) polynomial approximation of the inverse
    p = SVector{4}((x1, x2, x3, 1/df1, 1/df2, 1/df3))
    M = transpose(SMatrix{4, 4}(
        1.0,  f1, f1^2, f1^3, f1^4, f1^5,
        1.0,  f2, f2^2, f2^3, f2^4, f2^5,
        1.0,  f3, f3^2, f3^3, f3^4, f3^5,
        0.0, 1.0, 2*f1, 3*f1^2, 4*f1^3, 5*f1^4,
        0.0, 1.0, 2*f2, 3*f2^2, 4*f2^3, 5*f2^4,
        0.0, 1.0, 2*f3, 3*f3^2, 4*f3^3, 5*f3^4,
    ))
    # interpolate the inverse
    # f^(-1)(x) = b0 + b_i x^i, so
    # f^(-1)(0) = b0
    b0, _, _, _ = pinv(M)*p
    return b0
end