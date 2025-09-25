""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter)
end

""" 
Keplerian orbit evolution. Uses the universal variable formaulation from Danby, with initial guesses from Vallado, and
modification to the stumpff functions from Wisdam + Hernandez 2015
"""
function propagate(pos, vel, dt, gm; max_iter = 20)
    if all(vel .== 0) || all(pos .== 0) 
        throw("invalid orbit")
    end

    if dt == 0
        return (pos, vel)
    elseif dt < 0 # THE STUPID GUESS DOESNT WORK FOR BACKWARDS TIME
        posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf
    end

    # normalize so that r0 = 1, gm = 1
    # these are chosen to elliminate some computations
    # and so that the numbers we're working with are (hopefully) in a nice floating point regime
    DU   = norm(pos)
    TU   = sqrt(DU^3 / gm)
    pos0 = pos/DU
    vel0 = vel/DU*TU
    dt0  = dt/TU



    # println("recovered gm = $(DU^3 / TU^2)")

    # r0  = norm(pos0)
    # dr0 = dot(vel0, pos0) / r0
    # v02 = dot(vel0, vel0)
    # alpha = 2gm/r0 - v02
    dr0 = dot(vel0, pos0)
    a   = 2.0 - dot(vel0, vel0) # 1/semimajor axis! alpha = beta in these coordinates
    # println("a  = $(1/a*DU)")
    # println("a0 = $(1/(2/norm(pos)-norm(vel)^2/gm))")

    Hvec = cross(pos0, vel0)
    # Evec = cross(vel0, Hvec)/gm - pos0/r0
    Evec = cross(vel0, Hvec) - pos0
    e    = norm(Evec)
    # println("e = $e")

    
    # we have non-dimentionalized in a way that s = x, so we can mix equations from both Vallado and Danby
    x = if abs(a) < 1e-6 
        # near-parabolic
        # converted from Vallado.
        # need to study this
        # NOT a upper or lower bound, so still need to establish a bracket
        p = dot(Hvec, Hvec)
        s = 0.5acot(3dt*sqrt(1/p^3))
        w = atan(tan(s)^(1/3))
        2sqrt(p)*cot(2w)
    elseif a > 0
        # elliptic
        # since dx/dt = 1/r, the minimum radius yields the maximum rate of change.
        # so we can establish an upper bound on x
        dt0*a/(1 - e)
    elseif a < 0
        # hyperbolic
        # converted from Vallado.
        # same caveats as near-parabolic case
        sqrt(-1/a)*log(-2.0a*dt0 / (dr0 + sqrt(-1/a)*(1.0 - a)))
    else
        throw("ERROR 1/a = $a")
    end
    # println("initial x = $(x*sqrt(DU))")

    # check we have not overflowed our numbers precision (excessively long timespan)
    # @assert 0 < s < typemax(T)

    # set up the bracket
    # only handle the forward time case (due to the above check)
    x_br = 0.0
    y_br = -dt0

    br_step = x

    _, c1, c2, c3 = stumpff(a*x^2)
    y = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
    while sign(y) == sign(y_br)
        # take forward steps, maintaining the size of our interval
        # we could use some newton steps here, and would probably work better
        # but I'm lazy, and this is exceedingly likely to require at most 1 step
        x_br = x
        y_br = y
        x    += br_step
        _, c1, c2, c3 = stumpff(a*x^2)
        y    = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
    end
    # should add a safety check here

    # if sign(y) == sign(y_br)
    #     throw("bracket search failed => [$y_br, $y]")
    # end

    # br   = 0.0
    # br_y = universal_kepler(br, alpha, r0, dr0, gm) - dt
    # y    = universal_kepler(s, alpha, r0, dr0, gm) - dt
    # while sign(y) == sign(br_y)
    #     br   = s
    #     br_y = y
    #     s    = 2s
    #     y    = universal_kepler(s, alpha, r0, dr0, gm) - dt
    # end

    # we now have a bracket, and can find the root.
    # chain a couple of closures together to get the right args to the right places
    # TODO: implement a bracketed root finder using higher derivatives
    bracket = (x_br, x)
    # println("bracket = $bracket")
    _f1 = _x -> (_x, stumpff(a*_x^2)...)
    _f2 = ((_x, _, _c1, _c2, _c3),) -> _x*_c1 + dr0*_c2*_x^2 + _c3*_x^3
    # x = find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, A42())
    x = find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, Bisection())
    # x = try
    #     find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, Bisection())
    # catch err
    #     throw("bracket = $bracket")
    # end
    # println("final x = $(x*sqrt(DU))")
    # println("final z = $(a*x^2)")

    # _, c1, c2, c3 = stumpff(alpha*s^2)
    _, c1, c2, c3 = stumpff(a*x^2)
    # dtf = x*c1 + dr0*c2*x^2 + c3*x^3
    # println("initial dt = $(dt0)")
    # println("final dt   = $(dtf)")
    # println("dt error   = $((dtf-dt0))")

    # f = 1 - (gm/r0)*s^2*c2
    # g = dt - gm*s^3*c3
    # posf = f*pos0 + g*vel0

    # rf = norm(posf)
    # df = -(gm/(rf*r0))*s*c1
    # # dg = (r0/rf)*(c0 + dr0*c1)
    # dg = 1 - (gm/rf)*s^2*c2
    # velf = df*pos0 + dg*vel0

    f    = 1.0 - c2*x^2
    g    = dt0 - c3*x^3
    posf = f*pos0 + g*vel0
    # println("f = $(f)")
    # println("g = $(g*TU)")

    rf   = norm(posf)
    # df   = -c1*x/rf
    dg   = 1.0 - c2/rf*x^2
    df   = x*(a*c3*x^2 - 1.0)/rf
    velf = df*pos0 + dg*vel0
    # println("df = $(df)")
    # println("dg = $(dg)")

    # println("f*dg - df*g = $(f*dg - df*g)")

    # re-dimensionalize our inputs
    return posf*DU, velf*DU/TU
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

# dimensionalized initital guesses
    # s = if abs(alpha/gm) < parabolic_tol # (near) parabolic case
    #     # @debug "(near) parabolic orbit (alpha = $(alpha))"
    #     # solve the cubic solution to the parabolic case
    #     # a0 = -6dt/gm
    #     # a1 = 6r0/gm
    #     # a2 = 3r0*dr0/gm

    #     # q1 = a1/3 - (a2^2)/9
    #     # q2 = (a1*a2 - 3a0)/6 - (a2^3)/27

    #     # @assert q1^3 > q2^2

    #     # p1 = (q2 + sqrt(q1^3 + q2^2))^(1/3)
    #     # p2 = (q2 - sqrt(q1^3 + q2^2))^(1/3)

    #     # p1 + p2 - a2/3

    #     # from Vallado
    #     p = dot(Hvec, Hvec) / gm
    #     s = 0.5acot(3dt*sqrt(gm/p^3))
    #     w = atan(tan(s)^(1/3))
    #     2sqrt(p)*cot(2w)/sqrt(gm)
    # elseif alpha > 0.0 # elliptic
    #     # @debug "elliptic orbit (alpha = $(alpha))"
    #     dt*alpha/(gm*(1 - e))
    # else # hyperbolic
    #     # trying new hyperbolic guess 
    #     # @debug "hyperbolic orbit (alpha = $(alpha))"
    #     # k  = 1.8
    #     # dM = sqrt(-gm*alpha^3)*dt
    #     # CH = 1 - r0*alpha
    #     # SH = r0*dr0*sqrt(-alpha/gm)
    #     # dF = dM > 0 ? log((2dM + k*e)/(CH + SH)) : log((-2dM + k*e)/(CH - SH)) 
    #     # dF / sqrt(-alpha)
    #     a = alpha
    #     x0 = sign(dt)*sqrt(-1.0/a)*log(-2.0*gm*a*dt / (r0*dr0 + sign(dt)*sqrt(-gm*a)*(1.0 - r0*a)))
    #     x0/sqrt(gm)
    # end