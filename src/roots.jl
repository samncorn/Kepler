""" Chandrupatla's method, with an additional criteria taken from brent's method to guaruntee convergence
"""
function chandrupatla_brent(f, bracket)
    a, b = bracket
    ya = f(a)
    yb = f(b)
    @assert sign(ya) != sign(yb) "interval cannot be guarunteed to bracket a root"
    x  = a + (b - a)/2
    c  = a # contrapoint
    yc = ya
    println("a = $a, b = $b, x = $x")
    println("ya = $ya, yb = $yb")
    println()
    dx1 = Inf
    dx2 = Inf
    # dx3 = Inf
    i  = 1
    xp = x + eps(x)
    while x != xp
        y = f(x)
        if y == 0
            return x
        end
        println(i)
        println("a = $a, b = $b, x = $x")
        println("ya = $ya, yb = $yb, y = $y")
        println()
        if sign(y) == sign(yb)
            c = b
            b = a
            yc = yb
            yb = ya
        elseif sign(y) == sign(ya)
            # b becomes the contrapoint
            c  = a
            yc = ya
        else
            throw("bad value (ya = $ya, yb = $yb, y = $y) encountered")
        end
        a  = x
        ya = y

        xp  = x

        phi = (ya - yb) / (yc - yb)
        xi  = (a - b) / (c - b)
        println("ξ = $xi")
        # if dc <= da/2 && 
        # @assert a < b
        lim1 = phi^2
        lim2 = 1 - (1 - phi)^2
        lim1, lim2 = sort((lim1, lim2))
        if abs(c - b) < dx1/2 && lim1 < xi < lim2
        # if abs(c - b) < dx1/2 && xi > phi^2 && xi < 1 - (1 - phi)^2
            # inverse interpolate
            println("interpolate")
            k = (c - a) / (b - a)
            t = (ya/(ya - yb))*(yc/(yc - yb)) - k*(ya/(yc - ya))*(yb/(yb - yc))
            x = a + t*(b - a)
        else
            # bisect
            println("bisect")
            x = a + (b - a)/2
        end
        dx1 = dx2
        dx2 = abs(c - b)
        i += 1
    end
    return x
end

""" uses inverse quintic interpolation (3 points, 3 derivatives) to compute the next iteration. Falls back 
to inverse quadratic (3 points no derivatives), which falls back to a bisection. fall back criteria are both 
chandrupatla's criterion and Brents 
"""
function inverse_quintic(fdf, bracket)
    a, b = bracket
    ya, da = fdf(a)
    yb, db = fdf(b)
    @assert sign(ya) != sign(yb) "interval cannot be guarunteed to bracket a root"
    x  = a + (b - a)/2
    c  = a # contrapoint
    yc = ya
    dc = da
    dx1 = Inf
    dx2 = Inf
    i  = 1
    xp = x + eps(x)
    while x != xp
        y, dy = fdf(x)
        if y == 0
            return x
        end
        if sign(y) == sign(yb)
            c = b
            b = a
            yc = yb
            yb = ya
            dc = db
            db = da
        elseif sign(y) == sign(ya)
            # b becomes the contrapoint
            c  = a
            yc = ya
            dc = da
        else
            throw("bad value (ya = $ya, yb = $yb, y = $y) encountered")
        end
        a  = x
        ya = y
        da = dy

        xp = x

        phi = (ya - yb) / (yc - yb)
        xi  = (a - b) / (c - b)

        lim1 = phi^2
        lim2 = 1 - (1 - phi)^2
        lim1, lim2 = sort((lim1, lim2))
        if abs(c - b) < dx1/2 && lim1 < xi < lim2
            if da != 0 && db != 0 && dc != 0 && ya != yb && yb != yc && ya != yc # && phi^5 < xi < 1 - (1 - phi)^5
                # try our fancy pants inverse quintic approximation
                # x = inverse_quintic_step(a, b, c, ya, yb, yc, da, db, dc)
                x = inverse_quintic_step(c, b, a, yc, yb, ya, dc, db, da)
            else
                # fallback to inverse quadratic
                k = (c - a) / (b - a)
                t = (ya/(ya - yb))*(yc/(yc - yb)) - k*(ya/(yc - ya))*(yb/(yb - yc))
                x = a + t*(b - a)
            end
        else
            # bisect
            x = a + (b - a)/2
        end
        dx1 = dx2
        dx2 = abs(c - b)
        i += 1
        println(i)
    end
    return x
end

# overload for having separate function and derivative evaluations
function inverse_quintic(f, df, bracket)
    fdf = x -> (f(x), df(x))
    return inverse_quintic(fdf, bracket)
end

# function inverse_cubic_step(xn, xn1, yn, yn1, dn, dn1)

# end

function inverse_quintic_step(xn, xn1, xn2, yn, yn1, yn2, dn, dn1, dn2)
    q0 = yn/yn2 
    q1 = yn1/yn2
    a0 = q1^2*(q0*(3 + 3q1 - 5q0) - q1) / ((q0 - 1)*(q0 - q1))^3
    a1 = q0^2*(q1*(5q1 - 3q0 - 3) + q0) / ((q1 - 1)*(q0 - q1))^3
    a2 = (q0*q1)^2*(3q1 - q0*(q1 - 3) - 5) / ((q0 - 1)*(q1 -1))^3
    b0 = q0*q1^2 / ((q0 - 1)*(q0 - q1))^2
    b1 = q0^2*q1 / ((q0 - q1)*(q1 - 1))^2
    b2 = ((q0*q1) / ((q0 - 1)*(q1 - 1)))^2

    x4 = -yn2*(b0/dn + b1/dn1 + b2/dn2) - a0*xn - a1*xn1 - a2*xn2
    return x4
end

function bisect_step(xl, xh)
    return 0.5(xl + xh)
end

function lmm12_step(a, b, f1, f2, df1, df2)
    # construct polynomial approximation of the inverse
    # cubic in this case
    # p = SVector{4}((a, b, 1/df1, 1/df2))
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

    return h00*a + h10*(f2-f1)/df1 + h01*b + h11*(f2-f1)/df2
end

# function cubic_hermite_spline(x, a, b, ya, yb)


# function brent_lmm12_step(xa, del1, del2, a, b, f1, f2, df1, df2, previous_bisect; tol = 1e-15)
#     if previous_bisect && abs(del1) < tol
#         return 0.5(xa + a), true
#     end

#     if !previous_bisect && abs(del2) < tol
#         return 0.5(xa + a), true
#     end

#     s = lmm12_step(a, b, f1, f2, df1, df2)
#     if previous_bisect && 
#         return 0.5(xa + a), true
#     end

#     if !previous_bisect && !((3xa + a)/4 < s < a)
# end


function check_step(x, xl, xh; tol = 1e-14)
    return x == xl || x == xh || abs(xl - xh) < tol
end