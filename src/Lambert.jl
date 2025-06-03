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

    dm = if dm == Short
        1.0
    elseif dm == Long
        -1.0
    else
        throw("direction of motion must be short or long")
    end

    cosdv = dot(pos1, pos2) / (r1*r2)
    # sindv = dm*sqrt(1 - cosdv^2)

    A = dm*sqrt(r1*r2*(1 + cosdv))

    if A == 0
        return nothing
    end

    # we can set an inital bracket with the asymptotic boundaries for given n revs
    # 0 rev permits hyperbolic orbits, which have negative phi
    z1 = n == 0 ? -4*pi : (2*n*pi)^2
    z2 = (2*(n + 1)*pi)^2

    # but these initial guesses do not have defined values, much less derivatives
    # so the first couple of steps need to be bisections, until we establish an adequate bracket
    # also, since there are two solutions (either of high or low energy), we have to decide which
    # branch to follow

    # if de = High, dphi/dt < 0
    # if de = Low, dphi/dt > 0
    # so we can use this when handling the bracket

    # initialize rootfinding bracket
    # dt1, ddt1 = _universal_lambert_2(z1, A, r1, r2, gm)
    # dt2, ddt2 = _universal_lambert_2(z2, A, r1, r2, gm)
    # dt1, ddt1 = n == 0 ? _universal_lambert_2(z1, A, r1, r2, gm) : (NaN, NaN)
    dt1, ddt1 = (NaN, NaN)
    dt2, ddt2 = (NaN, NaN)

    if n == 0
        bracketed = false
        while !bracketed
            _, _, c2, c3 = stumpff(z1)
            y = r1 + r2 + A*(z1*c3 - 1)/sqrt(c2)

            if (A > 0 && y < 0)
                z1 *= 0.5
                _, _, c2, c3 = stumpff(z1)
                y = r1 + r2 + A*(z1*c3 - 1)/sqrt(c2)
            else
                dt1, ddt1 = _universal_lambert_2(z1, A, r1, r2, gm)
                if dt1 < dt
                    bracketed = true
                else
                    z1 *= 1.5
                end
            end
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

    println()
    println(i)
    println("$(z1), $(z2)")
    println("$(dt1 - dt), $(dt2 - dt)")

    # iterate until we establish a well-defined bracket
    while i < max_iter && err_abs > tol && err_rel > tol
        i   += 1
        dz2 = dz1
        dz1 = abs(zi - z)

        # err_abs = abs(dt1 - dt)
        # err_rel = err_abs

        dti, ddti = _universal_lambert_2(z, A, r1, r2, gm)
        if (dti - dt) <= tol
            # check if we found the wrong root
            if (de == High && ddti <= 0) || (de == Low && ddti >= 0)
                z1  = z
                z2  = z
                err_abs = 0.0
                err_rel = 0.0
                continue
            end
        end

        err_abs = abs(dti - dt)
        err_rel = err_abs

        # decide which side of the bracket to adjust
        if de == Low || n == 0
            if dti < dt || ddti < 0
                # make lower bound
                z1   = z
                dt1  = dti
                ddt1 = ddti
            else
                # make upper bound
                z2   = z
                dt2  = dti
                ddt2 = ddti
            end
        elseif de == High
            if dti < dt || ddti > 0
                # make upper bound
                z2   = z
                dt2  = dti
                ddt2 = ddti
            else
                # make lower bound
                z1   = z
                dt1  = dti
                ddt1 = ddti
            end
        end

        println()
        println("i = $(i)")
        # println("  z1 = $(z1) z2 = $(z2)")
        @printf "      z1 = % 10.7e     z2 = % 10.7e\n" z1 z2
        @printf "  dt1-dt = % 10.7e dt2-dt = % 10.7e\n" (dt1-dt) (dt2 - dt)
        @printf "      dz = % 10.7e\n" dz1

        # decide what kind of step to take
        zi = z
        if (isnan(dt1) || isnan(dt2)) || (dt1 > dt && dt2 > dt)
            @printf "       bisect\n"
            # still trying to properly bracket, try newton or bisection
            z = 0.5(z1 + z2)
            bs = true
            # if isnan(dt1) && isnan(dt2)
            #     println("  bisecting")
            #     # bisect
            #     z = 0.5(z1 + z2)
            #     bs = true
            # else
            #     println("  newton step")
            #     # newton step 
            #     if isnan(dt1)
            #         z = z2 - (dt2 - dt)/ddt2
            #     elseif isnan(dt2)
            #         z = z1 - (dt1 - dt)/ddt1
            #     end
            #     bs = false
            # end
        elseif dt1 < dt && dt2 < dt
            # something has gone horribly wrong, and both our endpoints are between the two roots
            throw("bracket between the two roots")
        else
            # have differentiable endpoints and a single root between them
            # so we can take a brent-newton step
            z, bs = brent_newton_step((z1, z2), (dt1-dt, dt2-dt), (ddt1, ddt2), dz1, dz2, bs)
            if bs 
                @printf "       bisect\n"
            else
                @printf "       interpolate\n"
            end
        end
        @printf "       z = % 10.7e\n" z
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

""" Computes both dt and dt/dphi
"""

function _universal_lambert_2(z, A, r1, r2, gm)
    _, _, c2, c3 = stumpff(z)
    y   = r1 + r2 + A*(z*c3 - 1)/sqrt(c2)
    s   = sqrt(y/(gm*c2))
    # dt  = gm*(c3 + A*sqrt(y))
    dt  = gm*c3*s^3 + A*sqrt(y/gm)
    dc2 = (1 - z*c3 - 2c2)/(2z)
    dc3 = (c2 - 3c3)/(2z)
    ddt = (dc3 - 3c3*dc2/(2c2))*gm*s^3 + (3c3*sqrt(y/gm)/c2 + A/(gm*s))*A/8
    return dt, ddt
end

# end