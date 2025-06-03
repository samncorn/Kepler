
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
    max_iters = ceil(Int, log2(tol/dx1))^2
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

# @enum root_result found_root bracket_error convergence_error

# """ f must return the appropriate amount of derivatives for next_term
# """
# function root_solve(next_term, f, x0, tol, max_iter)
#     i = 0
#     x = x0
#     err = Inf
#     while err > tol && i <= max_iter
#         dx = next_term(f, x)
#         x -= dx
#         err = abs(dx)
#         i += 1
#     end
#     return x, i
# end

# function _newton2(f, x)
#     y, dy = f(x)
#     dx = y/dy
#     return dx
# end

# function _newton3(f, x)
#     y, dy, ddy = f(x)
#     dx = y/dy
#     dx2 = y/(dy - dx*ddy/2)
#     return dx2
# end

# function _newton4(f, x)
#     y, dy, ddy, dddy = f(x)
#     dx = y/dy
#     dx2 = y/(dy - dx*ddy/2)
#     dx3 = y/(dy - dx2*ddy/2 + dx2^2 * dddy/6)
#     return dx3
# end

# function _laguerre(f, x, n)
#     y, dy, ddy = f(x)
#     return n*y / (dy + sign(dy)*sqrt(abs((n-1)^2 * (dy^2 - n*(n-1)*y*ddy))))
# end

# newton2(f, x0, tol, max_iter) = root_solve(_newton2, f, x0, tol, max_iter)
# newton3(f, x0, tol, max_iter) = root_solve(_newton3, f, x0, tol, max_iter)
# newton4(f, x0, tol, max_iter) = root_solve(_newton4, f, x0, tol, max_iter)
# # laguerre(f, x0, tol, max_iter; n = 5) = root_solve(x -> _laguerre(x[1], x[2], 5), f, x0, tol, max_iter)
# function laguerre(f, x0, tol, max_iter; n = 5)
#     i = 0
#     x = x0
#     err = Inf
#     while err > tol && i <= max_iter
#         dx = _laguerre(f, x, n)
#         x -= dx
#         err = abs(dx)
#         i += 1
#     end
#     return x, i
# end

# function _bracket_root(method, f, bracket, tol, max_iter)
#     # check the bracket, if no root return failure
#     if sign(f(bracket[1])) == sign(f(bracket[2]))
#         return bracket_err, bracket[1]
#     end

#     iter = 0
#     err  = abs(bracket[1] - bracket[2])

#     while err > tol && iter < max_iter
#         # cubic spline the inverse
#         x = 
    
#         err   = abs(bracket[1] - bracket[2])
#         iter += 1
#     end
#     return found_root, x
# end

# function _newton_bracket_2()

# end

# newton_bracket_2(f, bracket, tol, max_iter) = _bracket_root(_newton_bracket_2, f, bracket, tol, max_iter)