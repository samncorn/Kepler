# module Lambert

# using Roots
# using Stumpff

@enum DirectionOfEnergy High Low
@enum DirectionOfMotion Short Long

"""
analytically solve Lamberts problem. Algorithm from Vallado, 5th edition, pgs 504-505

If no valid solution exists, returns nothing
"""
function lambert_solve(pos1, pos2, dt, gm; de = High, dm = Short, n::Int = 0, max_iter = 1_000, tol = 1e-12)
    @assert n >= 0
    if dt < 0
        return lambert_solve(pos2, pos1, -dt, gm; de = de, dm = dm, n = n)
    end

    r1 = norm(pos1)
    r2 = norm(pos2)

    # dm = if dm == Short
    #     1.0
    # elseif dm == Long
    #     -1.0
    # else
    #     throw("direction of motion must be short or long")
    # end

    cosdv = dot(pos1, pos2) / (r1*r2)
    # sindv = dm*sqrt(1 - cosdv^2)

    # A = dm*s
    A = qrt(r1*r2*(1 + cosdv))
    A = dm == short ? A : -A

    if A == 0
        return nothing
    end

    # compute the value and first two derivatives of the lambert universal phi -> time function
    _dt(x)   = universal_lambert(x, A, r1, r2, gm)
    _ddt(x)  = Enzyme.gradient(ForwardWithPrimal, _dt, x)
    _dddt(x) = Enzyme.gradient(ForwardWithPrimal, _ddt, x)

    z1, z2 = lambert_find_bracket(A, r1, r2, gm, dt; n = 0, de = de, dm = dm)
    # we now have a bracket around a single root, although one end is asymptotic, so we need to bisect until we have a 
    # complete bracket
    # we also need to check if we accidentally found our root with the bracket

    # # dt1, ddt1 = _universal_lambert_2(z1, A, r1, r2, gm)
    # # dt2, ddt2 = _universal_lambert_2(z2, A, r1, r2, gm)
    # # dt1, ddt1 = n == 0 ? _universal_lambert_2(z1, A, r1, r2, gm) : (NaN, NaN)
    dt1, ddt1 = (NaN, NaN)
    dt2, ddt2 = (NaN, NaN)

    if n == 0
        derivs = _ddt(z1)
        dt1  = derivs.val 
        ddt1 = derivs.derivs[1]
    else
        if de == High
            derivs = _ddt(z2)
            dt2  = derivs.val 
            ddt2 = derivs.derivs[1]
        else
            derivs = _ddt(z1)
            dt1  = derivs.val 
            ddt1 = derivs.derivs[1]
        end
    end

    # use Brent-Newton to find the root
    i = 0
    err_abs = Inf
    err_rel = Inf
    bs  = true
    dz1 = abs(z1 - z2)
    dz2 = Inf

    # initialize the first step
    zi = z1
    z  = 0.5(z1 + z2)

    # println()
    # println("i = $(i)")
    # @printf "      z1 = % 10.7e     z2 = % 10.7e\n" z1 z2
    # @printf "  dt1-dt = % 10.7e dt2-dt = % 10.7e\n" (dt1-dt) (dt2 - dt)
    # @printf "      dz = % 10.7e\n" dz1

    # iterate until we establish a well-defined bracket
    while i < max_iter && err_abs > tol && err_rel > tol
        i   += 1
        # println()
        # println("i = $(i)")
        dz2  = dz1
        dz1 = abs(zi - z)

        derivs = _ddt(z)
        dti  = derivs.val
        ddti = derivs.derivs[1]

        err_abs = abs(dti - dt)
        err_rel = err_abs

        if abs(dti - dt) < tol
            z1   = z
            dt1  = dti
            ddt1 = ddti

            z2   = z
            dt2  = dti
            ddt2 = ddti

            err_abs = 0.0
            err_rel = 0.0
            # @printf "      z1 = % 10.7e     z2 = % 10.7e\n" z1 z2
            # @printf "  dt1-dt = % 10.7e dt2-dt = % 10.7e\n" (dt1-dt) (dt2 - dt)
            # @printf "      dz = % 10.7e\n" dz1
            continue
        end

        # decide which side of the bracket to adjust
        if (!isnan(dt1) && sign(dti - dt) == sign(dt1 - dt)) || (!isnan(dt2) && sign(dti - dt) != sign(dt2 - dt))
            # @printf "  left\n"
            z1   = z
            dt1  = dti
            ddt1 = ddti
        elseif !isnan(dt2) && sign(dti - dt) == sign(dt2 - dt) || (!isnan(dt1) && sign(dti - dt) != sign(dt1 - dt))
            # @printf "  right\n"
            z2   = z
            dt2  = dti
            ddt2 = ddti
        else
            throw("uncovered case")
        end

        # @printf "      z1 = % 10.7e     z2 = % 10.7e\n" z1 z2
        # @printf "  dt1-dt = % 10.7e dt2-dt = % 10.7e\n" (dt1-dt) (dt2 - dt)
        # @printf "      dz = % 10.7e\n" dz1

        # decide what kind of step to take
        zi = z
        if isnan(dt1) || isnan(dt2)
            # @printf "       bisect\n"
            # still trying to properly bracket, try newton or bisection
            z = 0.5(z1 + z2)
            bs = true
        elseif dt1 < dt && dt2 < dt
            # something has gone horribly wrong, and both our endpoints are between the two roots
            throw("bracket between the two roots")
        else
            # have differentiable endpoints and a single root between them
            # so we can take a brent-newton step
            z, bs = brent_newton_step((z1, z2), (dt1-dt, dt2-dt), (ddt1, ddt2), dz1, dz2, bs)
            if bs 
                # @printf "       bisect\n"
            else
                # @printf "       interpolate\n"
            end
        end
        # @printf "       z = % 10.7e\n" z
    end

    # compute fg coefficients
    _, _, c2, c3 = stumpff(z)
    y  = r1 + r2 + A*(z*c3 - 1)/sqrt(c2)
    f  = 1 - y/r1
    g  = A*sqrt(y/gm)
    dg = 1 - y/r2

    vel1 = (pos2 - f*pos1)/g
    vel2 = (dg*pos2 - pos1)/g

    return vel1, vel2
end

function lambert_find_bracket(A, r1, r2, gm, dt; n = 0, de = High, dm = Short, tol = 1e-12, max_iter = 100)
    @assert n >= 0
    z1 = n == 0 ? -4*pi : (2*n*pi)^2
    z2 = (2*(n + 1)*pi)^2

    _dt(x)   = universal_lambert(x, A, r1, r2, gm)
    _ddt(x)  = Enzyme.gradient(ForwardWithPrimal, _dt(x), x)
    _dddt(x) = Enzyme.gradient(ForwardWithPrimal, _ddt(x), x)

    # establish either an upper or lower defined (non-asymptotic) bound 
    if n == 0
        # zero-rev, find lower bound
        if dm == Short
            # in this case, dt(z) is udnefined for z < 0
            z1 = 0.0
        else
            # step until we have either found dti < dt, or dt(z1) is undefined, in which case 
            # we still have a bracket, we just have to bisect
            dz = 4pi
            z1 = 0.0
            bracketed = false
            while !bracketed
                _, _, c2, c3 = stumpff(z1)
                y = r1 + r2 + A*(z1*c3 - 1)/sqrt(c2)

                if (A > 0 && y < 0)
                    # exceeded domain, can return and bisect from there
                    bracketed = true
                else
                    dt1    = _dt(z1)
                    if dt1 <= dt
                        bracketed = true
                    else
                        z2 -= dz
                        z1 -= dz
                        # z1 *= 1.5
                    end
                end
            end
        end
    else
        # multi-rev, find local minimum (branch switch point)
        zi = z1
        z = 0.5(z1 + z2)
        zz1 = z1
        zz2 = z2

        dz2 = Inf
        dz1 = abs(z1 - z2)

        f1   = NaN
        df1  = NaN
        ddf1 = NaN

        f2   = NaN
        df2  = NaN
        ddf2 = NaN
        
        bs = false
        i  = 0

        derivs = _dddt(z)
        f0   = derivs.val.val
        df0  = derivs.val.deriv
        ddf0 = derivs.derivs[1].derivs[1]

        while i < max_iter && abs(df0) > tol
            i += 1

            dz2 = dz1
            dz1 = abs(z - zi)
            if df0 > 0
                # dz1  = (z - zz2)
                zz2  = z
                f2   = f0
                df2  = df0
                ddf2 = ddf0
            else
                # dz1  = (z - zz1)
                zz1  = z
                f1   = f0
                df1  = df0
                ddf1 = ddf0
            end

            if isnan(f1) || isnan(f2)
                # bisect
                z = 0.5(zz1 + zz2)
                bs = true
            else
                # newton-brent
                z, bs = brent_newton_step((zz1, zz2), (df1, df2), (ddf1, ddf2), dz1, dz2, bs)
            end

            derivs = _dddt(z)
            f0   = derivs.val.val
            df0  = derivs.val.deriv
            ddf0 = derivs.derivs[1].derivs[1]
        end

        if de == High
            z2 = z
        else
            z1 = z
        end
    end

    return z1, z2
end

function universal_lambert(z, A, r1, r2, gm)
    _, _, c2, c3 = stumpff(z)
    y   = r1 + r2 + A*(z*c3 - 1)/sqrt(c2)

    if (A > 0 && y < 0)
        return NaN
    end

    s   = sqrt(y/(gm*c2))

    dt  = gm*c3*s^3 + A*sqrt(y/gm)
    return dt
end

""" Computes both dt and dt/dphi
"""

# function _universal_lambert_2(z, A, r1, r2, gm)
#     _, _, c2, c3 = stumpff(z)
#     y   = r1 + r2 + A*(z*c3 - 1)/sqrt(c2)
#     s   = sqrt(y/(gm*c2))

#     dt  = gm*c3*s^3 + A*sqrt(y/gm)
#     dc2 = (1 - z*c3 - 2c2)/(2z)
#     dc3 = (c2 - 3c3)/(2z)
#     ddt = (dc3 - 3c3*dc2/(2c2))*gm*s^3 + (3c3*sqrt(y/gm)/c2 + A/(gm*s))*A/8
#     return dt, ddt
# end

# function _universal_lambert_3(z, A, r1, r2, gm)
#     _, _, c2, c3 = stumpff(z)
#     y   = r1 + r2 + A*(z*c3 - 1)/sqrt(c2)
#     s   = sqrt(y/(gm*c2))

#     dt   = gm*c3*s^3 + A*sqrt(y/gm)
#     dc2  = (1 - z*c3 - 2c2)/(2z)
#     dc3  = (c2 - 3c3)/(2z)
#     ddt  = (dc3 - 3c3*dc2/(2c2))*gm*s^3 + (3c3*sqrt(y/gm)/c2 + A/(gm*s))*A/8

#     q  = 
#     s1 = 
#     s2 = 
#     s3 = 
#     s4 = 
#     dddt = 
#     return dt, ddt, dddt
# end

# end