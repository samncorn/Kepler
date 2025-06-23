""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20, parabolic_tol = 1e-6, tol = 1e-12)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter, parabolic_tol = parabolic_tol, tol = tol)
end

""" 
Keplerian orbit evolution. Uses the universal variable formaulation from Danby, with initial guesses from Vallado, and
modification to the stumpff functions from Wisdam + Hernandez 2015
"""
function propagate(pos0, vel0, dt, gm; max_iter = 20, parabolic_tol = 1e-6, tol = 1e-12)
    if dt == 0
        return (pos0, vel0)
    elseif dt < 0 # THE STUPID GUESS DOESNT WORK FOR BACKWARDS TIME
        posf, velf = solve(pos0, -vel0, -dt, gm; max_iter = max_iter, parabolic_tol = parabolic_tol)
        return posf, -velf
    end

    r0  = norm(pos0)
    dr0 = dot(vel0, pos0) / r0
    v02 = dot(vel0, vel0)
    alpha = 2gm/r0 - v02

    Hvec = cross(pos0, vel0)
    Evec = cross(vel0, Hvec)/gm - pos0/r0
    e = norm(Evec)

    # initial guesses from Danby.

    s = if abs(alpha/gm) < parabolic_tol # (near) parabolic case
        # @debug "(near) parabolic orbit (alpha = $(alpha))"
        # solve the cubic solution to the parabolic case
        # a0 = -6dt/gm
        # a1 = 6r0/gm
        # a2 = 3r0*dr0/gm

        # q1 = a1/3 - (a2^2)/9
        # q2 = (a1*a2 - 3a0)/6 - (a2^3)/27

        # @assert q1^3 > q2^2

        # p1 = (q2 + sqrt(q1^3 + q2^2))^(1/3)
        # p2 = (q2 - sqrt(q1^3 + q2^2))^(1/3)

        # p1 + p2 - a2/3

        # from Vallado
        p = dot(Hvec, Hvec) / gm
        s = 0.5acot(3dt*sqrt(gm/p^3))
        w = atan(tan(s)^(1/3))
        2sqrt(p)*cot(2w)/sqrt(gm)
    elseif alpha > 0.0 # elliptic
        # @debug "elliptic orbit (alpha = $(alpha))"
        dt*alpha/(gm*(1 - e))
    else # hyperbolic
        # trying new hyperbolic guess 
        # @debug "hyperbolic orbit (alpha = $(alpha))"
        k  = 1.8
        dM = sqrt(-gm*alpha^3)*dt
        CH = 1 - r0*alpha
        SH = r0*dr0*sqrt(-alpha/gm)
        dF = dM > 0 ? log((2dM + k*e)/(CH + SH)) : log((-2dM + k*e)/(CH - SH)) 
        dF / sqrt(-alpha)
    end

    # check we have not overflowed our numbers precision (excessively long timespan)
    # @assert 0 < s < typemax(T)

    # set up the bracket
    # only handle the forward time case (due to the above check)
    br   = 0.0
    br_y = universal_kepler(br, alpha, r0, dr0, gm) - dt
    y    = universal_kepler(s, alpha, r0, dr0, gm) - dt
    while sign(y) == sign(br_y)
        br   = s
        br_y = y
        s    = 2s
        y    = universal_kepler(s, alpha, r0, dr0, gm) - dt
    end
    # @assert 0 < s < typemax(T)

    bracket = (br, s)
    s = find_zero(x -> universal_kepler(x, alpha, r0, dr0, gm) - dt, bracket, A42())

    # # chain some closures together
    # _f = x -> universal_kepler2(x, alpha, r0, dr0, gm)
    # _g = x -> (x[1] - dt, x[2])
    # _h = x -> _g(_f(x))
    # try
    #     s  = brent_newton(_h, bracket; max_iters = max_iter, tol = tol)
    # catch e
    #     throw("failed with error $(e) for α = $(alpha)")
    # end

    # @debug "s: $(s)"

    _, c1, c2, c3 = stumpff(alpha*s^2)

    f = 1 - (gm/r0)*s^2*c2
    g = dt - gm*s^3*c3
    posf = f*pos0 + g*vel0

    rf = norm(posf)
    df = -(gm/(rf*r0))*s*c1
    # dg = (r0/rf)*(c0 + dr0*c1)
    dg = 1 - (gm/rf)*s^2*c2
    velf = df*pos0 + dg*vel0

    return posf, velf
end

function universal_kepler(s, alpha, r0, dr0, gm)
    _, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + gm*s^3*c3
    return dt
end

function universal_kepler2(s, alpha, r0, dr0, gm)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + gm*s^3*c3
    r  = r0*c0   + r0*dr0*s*c1   + gm*s^2*c2
    return dt, r
end

function universal_kepler3(s, alpha, r0, dr0, gm)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + gm*s^3*c3
    r  = r0*c0   + r0*dr0*s*c1   + gm*s^2*c2
    dr = (-r0*alpha + gm)*s*c1 + r0*dr0*c0
    # dr = r0*dr0*c0 + gm*(1 - alpha*r0)*s*c1
    return dt, r, dr
end

function universal_kepler4(s, alpha, r0, dr0, gm)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt  = r0*s*c1 + r0*dr0*s^2*c2 + gm*s^3*c3
    r   = r0*c0   + r0*dr0*s*c1   + gm*s^2*c2
    dr  = (-r0*alpha + gm)*s*c1 + r0*dr0*c0
    ddr = (-r0*alpha + gm)*c0   - r0*dr0*alpha*s*c1
    # ddr = k*((1 - alpha*r0)*c0 - alpha*r0*dr0*s*c1)
    return dt, r, dr, ddr
end