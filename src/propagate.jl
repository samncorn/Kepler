""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter)
end

"Universal kepler solver."
function propagate(pos, vel, dt, gm; max_iter::Int = 20, tol = 1e-15)
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

    # better initial guesses
    # xh = dt/r0
    xh = if abs(gm*b) < 1e-6
        # parabolic (Vallado)
        h = cross(pos, vel)
        p = dot(h, h)/gm
        s = acot(3*sqrt(gm/p^3)*dt)/2
        w = atan(cbrt(tan(s)))
        sqrt(p)*2*cot(2w)/sqrt(gm)
    elseif gm*b < 0
        # hyperbolic (Vallado)
        a = gm/b
        abs(sqrt(-a)*log(-2gm*dt/(a*(s0+sqrt(-gm*a)*(1 - r0/a))))/sqrt(gm))
    elseif gm*b > 0
        # elliptic
        dt/r0
    end

    dth, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh = dth - dt

    # establish the bracket
    i = 0
    while i < 1000 && (sign(yl) == sign(yh) || isinf(dth) || isnan(dth))
        i += 1
        if isinf(dth) || isnan(dth) #|| isnan(rh)
            xh = (xl + xh)/2
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        elseif sign(yl) == sign(yh) && !isnan(rh)
            xl  = xh
            # xh += (dt - dth)/rh
            xh *= 2
            yl  = yh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh  = dth - dt
        end
    end

    # we can now guaruntee a solution
    # TODO: utilize derivative based methods 
    # x = 
    while abs(xh - xl) > tol
        x = 0.5(xl + xh)
        # x = 
        y = universal_kepler(x, b, r0, s0, gm) - dt
        if sign(y) == sign(yl)
            xl = x
            yl = y
        elseif sign(y) == sign(yh)
            xh = x
            yh = y
        elseif sign(y) == 0.0
            break
        end
    end
    x = 0.5(xl + xh)

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

function propagate_with_partials(pos, vel, dt, gm; max_iter = 20, tol = 1e-15)
    if dt == 0
        return pos, vel, I3, I3, I3, I3
    end

    if dt < 0
        # posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        posf, velf, dxdx, dxdv, dvdx, dvdv = propagate_with_partials(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf, dxdx, -dxdv, -dvdx, dvdv
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

    # better initial guesses
    # xh = dt/r0
    xh = if abs(gm*b) < 1e-6
        # parabolic (Vallado)
        h = cross(pos, vel)
        p = dot(h, h)/gm
        s = acot(3*sqrt(gm/p^3)*dt)/2
        w = atan(cbrt(tan(s)))
        sqrt(p)*2*cot(2w)/sqrt(gm)
    elseif gm*b < 0
        # hyperbolic (Vallado)
        a = gm/b
        abs(sqrt(-a)*log(-2gm*dt/(a*(s0+sqrt(-gm*a)*(1 - r0/a))))/sqrt(gm))
    elseif gm*b > 0
        # elliptic
        dt/r0
    end

    dth, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh = dth - dt

    i = 0
    while i < 1000 && (sign(yl) == sign(yh) || isinf(dth) || isnan(dth))
        i += 1
        if isinf(dth) || isnan(dth) #|| isnan(rh)
            xh = (xl + xh)/2
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        elseif sign(yl) == sign(yh) && !isnan(rh)
            xl  = xh
            # xh += (dt - dth)/rh
            xh *= 2
            yl  = yh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh  = dth - dt
        end
    end

    # we can now guaruntee a solution
    # TODO: utilize derivative based methods 
    while abs(xh - xl) > tol
        x = 0.5(xl + xh)
        y = universal_kepler(x, b, r0, s0, gm) - dt
        if sign(y) == sign(yl)
            xl = x
            yl = y
        elseif sign(y) == sign(yh)
            xh = x
            yh = y
        elseif sign(y) == 0.0
            break
        end
    end
    x = 0.5(xl + xh)

    # compute f and g functions
    _, c1, c2, c3, c4, c5 = stumpff5(b*x^2)
    f    = 1 - (gm/r0)*(x^2)*c2
    g    = dt - gm*(x^3)*c3
    posf = f*pos + g*vel
    rf   = norm(posf)
    df   = -(gm/(rf*r0))*x*c1
    dg   = 1 - (gm/rf)*(x^2)*c2
    velf = df*pos + dg*vel

    # compute the partials (Battin)
    C = gm*(x^2)*(gm*(x^3)*(3c5 - c4) - dt*c2)
    delr = posf - pos
    delv = velf - vel

    # try the Der 1998 formulation
    v0 = norm(vel)
    a  = b/gm
    U1 = sqrt(gm)*x*c1
    U2 = gm*(x^2)*c2
    U3 = sqrt(gm)*(x^3)*gm*c3
    M1 = pos*transpose(pos)/r0^2
    M2 = pos*transpose(vel)/(r0*v0)
    M3 = vel*transpose(pos)/(r0*v0)
    M4 = vel*transpose(vel)/v0^2

    k11  = abs(a) < 1e-6 ? 0.0 : (1/(a*rf*r0^2))*(3U1*U3 + (a*r0 - 2)*U2^2) + (U1^2)/rf + U2/r0
    k12  = v0*U1*U2/(rf*sqrt(gm))
    k13  = abs(a) < 1e-6 ? 0.0 : v0/(a*rf*sqrt(gm)*r0^2)*(r0*U1*U2 + 2*s0/sqrt(gm)*U2^2 + 3*U2*U3 - 3rf*U3 + a*(r0^2)*U1*U2)
    k14  = (v0^2)*(U2^2)/(rf*gm)
    dxdx = f*I3 + k11*M1 + k12*M2 + k13*M3 + k14*M4

    k21  = r0*U1*U2/(rf*sqrt(gm))
    k22  = abs(a) < 1e-6 ? 0.0 : (v0/(a*rf*gm))*(3U1*U3 + (a*r0 - 2)*U2^2)
    k23  = r0*v0*(U2^2)/(rf*gm)
    k24  = abs(a) < 1e-6 ? 0.0 : ((v0^2)/(a*rf*gm*sqrt(gm)))*(r0*U1*U2 + 2s0/sqrt(gm)*U2^2 + 3U2*U3 - 3rf*U3)
    dxdv = g*I3 + k21*M1 + k22*M2 + k23*M3 + k24*M4

    # dxdx = f*I3 + (rf/gm)*delv*transpose(delv) + (1/r0^3)*(r0*(1-f)*posf*transpose(pos) + C*velf*transpose(pos))
    # dxdv = g*I3 + (r0/gm)*(1-f)*(delr*transpose(vel) - delv*transpose(pos)) + (C/gm)*velf*transpose(vel)
    dvdx = (
        - (1/r0^2)*delv*transpose(pos) 
        - (1/rf^2)*posf*transpose(delv) 
        + df*(I3 - (1/rf^2)*posf*transpose(posf) + (1/(rf*gm))*(posf*transpose(velf) - velf*transpose(posf))*posf*transpose(delv))
        - gm*C/(rf*r0)^3*posf*transpose(pos)
        )
    dvdv = (r0/gm)*delv*transpose(delv) + (1/rf^3)*(r0*(1-f)*posf*transpose(pos) - C*posf*transpose(vel)) + dg*I3

    return posf, velf, dxdx, dxdv, dvdx, dvdv
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

function stumpff5(z)
    if z > 1e-3
        c0 = cos(sqrt(z))
        c1 = sin(sqrt(z))/sqrt(z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        c4 = (1/2 - c2)/z
        c5 = (1/6 - c3)/z
        return c0, c1, c2, c3, c4, c5
    elseif z < -1e-3
        c0 = cosh(sqrt(-z))
        c1 = sinh(sqrt(-z))/sqrt(-z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        c4 = (1/2 - c2)/z
        c5 = (1/6 - c3)/z
        return c0, c1, c2, c3, c4, c5
    else
        # z very small. evaluate the series to 6 terms 
        # error is O(x^7) < 1e-21, well within floating point tolerances
        # c2 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/182)/132)/90)/56)/30)/12)/2
        # c3 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/210)/156)/110)/72)/42)/20)/6
        c4 = stumpff_continued(4, z)
        c5 = stumpff_continued(5, z)
        c2 = 1/2 - z*c4
        c3 = 1/6 - z*c5
        c0 = 1 - z*c2
        c1 = 1 - z*c3
        
        return c0, c1, c2, c3, c4, c5
    end
end

""" evaluate the ith stumpff function as a continued fraction. Uses a fixed number of terms, so should be restricted to the regime near 0 (z < 1e-3 for double precision).
"""
function stumpff_continued(i, z)
    n = 6
    c = one(z)
    for k in n:-1:1
        # c -= z/((i + 2k)*(i + k))
        c = 1 - c*z/((i + 2k)*(i + k))
    end
    return c*factorial(i)
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