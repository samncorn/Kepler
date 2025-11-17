""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter)
end

"Universal kepler solver. normalizes to canonical units"
function propagate(pos, vel, dt, gm; max_iter = 20)
    if dt == 0
        return pos, vel
    end

    if dt < 0
        posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf
    end

    # constants
    r0 = norm(pos)
    s0 = dot(pos, vel)
    b  = 2gm/r0 - dot(vel, vel)

    # try to bracket the root
    # dt/dx = r >= 0, and monotonic, so we can step until we find a bracket
    # bracket = (0.0, dt/r)
    xl = 0.0
    yl = -dt

    xh = dt/r0
    dth, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh = dth - dt

    while sign(yl) == sign(yh) || isinf(dth) || isnan(dth)
        if sign(yl) == sign(yh)
            xl  = xh
            xh += (dt - dth)/rh
            yl  = yh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh  = dth - dt
        elseif isinf(dth) || isnan(dth)
            xh = (xl + xh)/2
            dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        end
    end

    x = try
        find_zero(_x -> universal_kepler(_x, b, r0, s0, gm) - dt, (xl, xh), A42())
    catch _
        println("orbit failed, dumping initial conditions...")
        println("pos = $pos")
        println("vel = $vel")
        println("dt = $dt")
        println("gm = $gm")
        throw("bad root find")
    end

    # compute f and g functions
    _, c1, c2, c3 = stumpff(b*x^2)
    f    = 1 - (gm/r0)*(x^2)*c2
    g    = dt - gm*(x^3)*c3
    posf = f*pos + g*vel
    rf   = norm(posf)
    df   = -(gm/(rf*r0))*x*c1
    dg   = 1 - (gm/rf)*(x^2)*c2
    velf = df*pos + dg*vel

    return posf, velf
end

function propagate_with_partials(pos, vel, dt, gm; max_iter = 20)
    if dt == 0
        return pos, vel
    end

    if dt < 0
        posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf
    end

    # constants
    r0 = norm(pos)
    s0 = dot(pos, vel)
    b  = 2gm/r0 - dot(vel, vel)

    # try to bracket the root
    # dt/dx = r >= 0, and monotonic, so we can step until we find a bracket
    # bracket = (0.0, dt/r)
    xl = 0.0
    yl = -dt

    xh = dt/r0
    dth, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh = dth - dt

    while sign(yl) == sign(yh) || isinf(dth)
        if sign(yl) == sign(yh)
            xl  = xh
            xh += (dt - dth)/rh
            yl  = yh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh  = dth - dt
        elseif isinf(dth)
            xh = (xl + xh)/2
            dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        end
    end
  
    x = try
        find_zero(_x -> universal_kepler(_x, b, r0, s0, gm) - dt, (xl, xh), A42())
    catch _
        println("orbit failed, dumping initial conditions...")
        println("pos = $pos")
        println("vel = $vel")
        println("dt = $dt")
        println("gm = $gm")
        throw("bad root find")
    end

    # compute f and g functions
    _, c1, c2, c3 = stumpff(b*x^2)
    f    = 1 - (gm/r0)*(x^2)*c2
    g    = dt - gm*(x^3)*c3
    posf = f*pos + g*vel
    rf   = norm(posf)
    df   = -(gm/(rf*r0))*x*c1
    dg   = 1 - (gm/rf)*(x^2)*c2
    velf = df*pos + dg*vel

    # partials (Goodyear 1965)
    U = x^2 * (c2*dt + gm*x^5*(b*x*c4 - 3c5))
    W = g*c2*y^2 + 3*gm*U 
    # intermediate matrix
    M = SMatrix{3, 3}(
        (U/(rf*r0) + 1/r0^2 + 1/rf^2)df - W*gm^2/(rf*r0)^3,
        df*y*c1/rf + (dg - 1)/rf^2,
        (dg - 1)*y*c1/rf - W*gm/rf^3,
        -df*y*c1/r0 - (f - 1)/r0^2,
        -df*c2*y^2,
        -(dg - 1)*c2*y^2,
        (f - 1)*c1*y/r0 - W*gm/r0^3,
        (f - 1)*c2*y^2,
        g*c2*y^2 - W
    )

    # sub matrices
    RV  = hcat(posf, velf)
    RV0 = transpose(hcat(pos, vel))
    M1  = SMatrix{2, 2}(M[2, 1], M[2, 2], M[3, 1], M[3, 2])
    M2  = SMatrix{2, 2}(M[2, 2], M[2, 3], M[3, 2], M[3, 3])
    M3  = SMatrix{2, 2}(M[1, 1], M[1, 2], M[2, 1], M[2, 2])
    M4  = SMatrix{2, 2}(M[1, 2], M[1, 3], M[2, 2], M[2, 3])

    PHI11 = f*I3 + RV*M1*RV0
    PHI12 = g*I3 + RV*M2*RV0
    PHI21 = df*I3 - RV*M3*RV0
    PHI22 = dg*I3 - RV*M4*RV0

    # final 6x6
    # PHI = vcat(hcat(PHI11, PHI12), hcat(PHI21, PHI22))

    return posf, velf, PHI11, PHI12, PHI21, PHI22
end

# function universal_kepler(y, l, k1, k2, k3)
#     _, c1, c2, c3 = stumpff(l*y^2)
#     L = y*(k1*c1 + y*(k2*c2 + y*k3*c3))
#     return L
# end

function universal_kepler(x, b, r0, s0, gm)
    z  = b*x^2
    _, c1, c2, c3 = stumpff(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    return dt
end

function universal_kepler2(x, b, r0, s0, gm)
    z  = b*x^2
    c0, c1, c2, c3 = stumpff(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    r  = r0*c0 + x*(s0*c1 + x*gm*c2)
    return dt, r
end

function stumpff(z)
    if z > 1e-3
        c0 = cos(sqrt(z))
        c1 = sin(sqrt(z))/sqrt(z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        return c0, c1, c2, c3
    elseif z < -1e-3
        c0 = cosh(sqrt(-z))
        c1 = sinh(sqrt(-z))/sqrt(-z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        return c0, c1, c2, c3
    else
        # z very small. evaluate the series to 6 terms 
        # error is O(x^7) < 1e-21, well within floating point tolerances
        c2 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/182)/132)/90)/56)/30)/12)/2
        c3 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/210)/156)/110)/72)/42)/20)/6
        c0 = 1 - z*c2
        c1 = 1 - z*c3
        return c0, c1, c2, c3
    end
end
""" 
[DEPRECATED]
Compute the position/velcoity along a Keplerian orbit after given timespan. Uses the universal variable formaulation
from Danby, with initial guesses from Vallado, and modification to the stumpff functions from Wisdam + Hernandez 2015.

normalizes the orbit to r0 = 1 and gm = 1, so there will likely be slight differences with other implementations from 
floating point operations. This conveniently leads to s (Danby) equal to x (Vallado), and alpha (Vallado) equal to 
beta (Wisdom, alpha in Danby). Total number of operations is reduced (though not likely impactful to most users)

Uses a bracketing method to guaruntee convergence IF the upper bound guess is well behaved (the lower bound is 
guarunteed). It is possible for the initial guess to evaluate to infinity if the timespan is too long, or if an orbit 
is near enough to rectilinear that numerical issues lead to e = 1. While dropping something into the sun is perfectly
well defined up to impact time, the singularity in the equations of motion causes a problem afterwards. More 
importantly it invalidates the initial guess. There is an argument to be made that solve up to the point fof numeric 
problems is desirable, (since we solve long time spans to the point of failure), and inward falling orbits are 
physically reasonable, so future revisions may include such a solution. Until I find that solution, no rectilinear 
orbits. 
"""
function _propagate(pos, vel, dt, gm; max_iter = 20)
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

    dr0 = dot(vel0, pos0)
    a   = 2.0 - dot(vel0, vel0) # 1/semimajor axis! alpha = beta in these coordinates

    Hvec = cross(pos0, vel0)
    Evec = cross(vel0, Hvec) - pos0
    e    = norm(Evec)
    
    # we have non-dimentionalized in a way that s = x, so we can mix equations from both Vallado and Danby
    x = if abs(a) < 1e-6 
        # near-parabolic
        # converted from Vallado.
        # need to study this
        # NOT an upper bound, so still need to establish a bracket
        p = dot(Hvec, Hvec)
        s = 0.5acot(3dt*sqrt(1/p^3))
        w = atan(tan(s)^(1/3))
        2sqrt(p)*cot(2w)
    else
        x = dt0*a/(1 - e)
    end

    if x == Inf
        throw("ERROR initial guess Inf. Likely causes are long timespans or near-rectilinear orbits. a = $(1/a*DU) e = $e q = $(a*(1-e)) dt = $dt pos = $pos vel = $vel gm = $gm")
    end

    # TODO: Add a check for near-rectilinear orbits. devise a solution

    # set up the bracket
    # only handle the forward time case (due to the above check)
    x_br = 0.0
    y_br = -dt0

    # _, c1, c2, c3 = stumpff(a*x^2)
    # y = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
    _, g1, g2, g3 = stumpff(a, x)
    y = g1 + dr0*g2 + g3 - dt0
    
    while sign(y) == sign(y_br)
        # take forward steps, maintaining the size of our interval
        # we could use some newton steps here, and would probably work better
        # but I'm lazy, and this is exceedingly likely to require at most 1 step
        x_br = x
        y_br = y
        x    *= 2
        # _, c1, c2, c3 = stumpff(a*x^2)
        # y    = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
        _, g1, g2, g3 = stumpff(a, x)
        y = g1 + dr0*g2 + g3 - dt0
    end

    # we now have a bracket, and can find the root.
    # chain a couple of closures together to get the right args to the right places
    # TODO: implement a bracketed root finder using higher derivatives
    bracket = (x_br, x)
    # println("bracket = $bracket")
    # _f1 = _x -> (_x, stumpff(a*_x^2)...)
    # _f2 = ((_x, _, _c1, _c2, _c3),) -> _x*_c1 + dr0*_c2*_x^2 + _c3*_x^3
    _f1 = _x -> stumpff(a, _x)
    _f2 = ((_, _g1, _g2, _g3),) -> _g1 + dr0*_g2 + _g3
    x = find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, A42())

    # _, c1, c2, c3 = stumpff(a*x^2)
    g0, g1, g2, _ = stumpff(a, x)
    f    = 1.0 - g2
    g    = g1 + dr0*g2
    posf = f*pos0 + g*vel0

    rf = norm(posf)
    df = -g1/rf
    dg = (g0 + dr0*g1)/rf
    # dg   = 1.0 - g2/rf
    # df   = x*(a*c3*x^2 - 1.0)/rf
    velf = df*pos0 + dg*vel0

    # re-dimensionalize our inputs
    return posf*DU, velf*DU/TU
end