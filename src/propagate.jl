""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter)
end

"Universal kepler solver."
function propagate(pos, vel, dt, gm; max_iter::Int = 100_000)
    if dt == 0
        return pos, vel
    end

    if dt < 0
        posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf
    end

    DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    # constants
    # r0 = norm(pos)
    # s0 = dot(vel, pos)
    # b  = 2gm/r0 - dot(vel, vel)
    b = 2 - dot(vel, vel)

    # x = kepler_guess(pos, vel, dt, gm)
    # x = solve_kepler_universal_A42(pos, vel, gm, dt)
    x = solve_kepler_universal_A42_canonical(pos, vel, dt)
    # x = solve_kepler_universal_lmm12(pos, vel, gm, dt; max_iter = max_iter)

    # compute f and g functions
    _, c1, c2, c3 = stumpff(b*x^2)
    # _, c1, c2, c3, _, _ = stumpff_fold(b*x^2)

    # more numerically precise? According to Rein + Tamayo WHFast paper
    # f    = -(gm/r0)*(x^2)*c2
    # g    = dt - gm*(x^3)*c3
    # posf = f*pos + g*vel + pos
    f = -c2*x^2
    g = dt - c3*x^3
    posf = f*pos + g*vel + pos
    
    # rf   = norm(posf)
    # df   = -(gm/(rf*r0))*x*c1
    # dg   = -(gm/rf)*(x^2)*c2
    # velf = df*pos + dg*vel + vel
    rf = norm(posf)
    df = -x*c1/rf
    dg = -(c2/rf)*x^2
    velf = df*pos + dg*vel + vel

    return posf*DU, velf*DU/TU
end

function propagate_with_partials(pos, vel, dt, gm; max_iter = 100_000)
    if dt == 0
        return pos, vel, I3, I3, I3, I3
    end

    if dt < 0
        # posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        posf, velf, dxdx, dxdv, dvdx, dvdv = propagate_with_partials(pos, -vel, -dt, gm; max_iter = max_iter)
        return posf, -velf, dxdx, -dxdv, -dvdx, dvdv
    end

        DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    # constants
    # r0 = norm(pos)
    # s0 = dot(vel, pos)
    # b  = 2gm/r0 - dot(vel, vel)
    b = 2 - dot(vel, vel)

    # x = kepler_guess(pos, vel, dt, gm)
    # x = solve_kepler_universal_A42(pos, vel, gm, dt)
    x = solve_kepler_universal_A42_canonical(pos, vel, dt)
    # x = solve_kepler_universal_lmm12(pos, vel, gm, dt; max_iter = max_iter)

    # compute f and g functions
    # _, c1, c2, c3 = stumpff(b*x^2)
    _, c1, c2, c3, c4, c5 = stumpff5(b*x^2)

    # more numerically precise? According to Rein + Tamayo WHFast paper
    # f    = -(gm/r0)*(x^2)*c2
    # g    = dt - gm*(x^3)*c3
    # posf = f*pos + g*vel + pos
    f = -c2*x^2
    g = dt - c3*x^3
    posf = f*pos + g*vel + pos
    
    # rf   = norm(posf)
    # df   = -(gm/(rf*r0))*x*c1
    # dg   = -(gm/rf)*(x^2)*c2
    # velf = df*pos + dg*vel + vel
    rf = norm(posf)
    df = -x*c1/rf
    dg = -(c2/rf)*x^2
    velf = df*pos + dg*vel + vel

    # compute the partials (Battin)
    C = (x^2)*((x^3)*(3c5 - c4) - dt*c2)
    delr = posf - pos
    delv = velf - vel

    # dxdx = f*I3 + (rf/gm)*delv*transpose(delv) + (1/r0^3)*(r0*(1-f)*posf*transpose(pos) + C*velf*transpose(pos))
    # units of DU/DU
    dxdx = (1 + f)*I3 + rf*delv*transpose(delv) + (-f)*posf*transpose(pos) + C*velf*transpose(pos)

    # dxdv = g*I3 + (r0/gm)*(1-f)*(delr*transpose(vel) - delv*transpose(pos)) + (C/gm)*velf*transpose(vel)
    # units of DU/VU = DU/(DU/TU) = TU
    dxdv = g*I3 + (-f)*(delr*transpose(vel) - delv*transpose(pos)) + C*velf*transpose(vel)

    # r0 = 1.0
    # dvdx = (
    #     - (1/r0^2)*delv*transpose(pos) 
    #     - (1/rf^2)*posf*transpose(delv) 
    #     + df*(I3 - (1/rf^2)*posf*transpose(posf) + (1/(rf*gm))*(posf*transpose(velf) - velf*transpose(posf))*posf*transpose(delv))
    #     - gm*C/(rf*r0)^3*posf*transpose(pos)
    #     )
    # Units of 1/TU
    dvdx = (
        - delv*transpose(pos) 
        - (1/rf^2)*posf*transpose(delv) 
        + df*(I3 - (1/rf^2)*posf*transpose(posf) + (1/rf)*(posf*transpose(velf) - velf*transpose(posf))*posf*transpose(delv))
        - C/(rf)^3*posf*transpose(pos)
        )

    # Units of VU/VU = 1
    # dvdv = (r0/gm)*delv*transpose(delv) + (1/rf^3)*(r0*(1-f)*posf*transpose(pos) - C*posf*transpose(vel)) + dg*I3
    dvdv = delv*transpose(delv) + (1/rf^3)*((-f)*posf*transpose(pos) - C*posf*transpose(vel)) + (1 + dg)*I3

    return posf*DU, velf*DU/TU, dxdx, dxdv*TU, dvdx/TU, dvdv
end

# function comp_sum!(init, v)
#     c = zero(init)

# end

function solve_kepler_universal_A42_canonical(pos, vel, dt)
    # constants
    s0 = dot(vel, pos)
    b  = 2 - dot(vel, vel)

    # better initial guesses
    xh = kepler_guess_canonical(pos, vel, dt)
    # xh = dt

    yh, rh = universal_kepler2_canonical(xh, b, s0)
    yh -= dt

    xl = zero(xh)
    # xl = 0.0
    yl = -dt
    rl = one(rh)

    # establish the bracket
    i = 0
    while i < 1000 && (sign(yl) == sign(yh) || isinf(yh) || isnan(yh))
        i += 1
        if isinf(yh) || isnan(yh) #|| isnan(rh)
            xh = (xl + xh)/2
            yh, rh = universal_kepler2_canonical(xh, b, s0)
            yh -= dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        elseif sign(yl) == sign(yh) && !isnan(rh)
            xl  = xh
            # xh += (dt - yh)/rh
            xh *= 2
            yl  = yh
            rl  = rh
            yh, rh = universal_kepler2_canonical(xh, b, s0)
            yh -= dt
        end
    end

    return find_zero(_x -> universal_kepler_canonical(_x, b, s0) - dt, (xl, xh), A42())
end

function solve_kepler_universal_A42(pos, vel, gm, dt)
    # constants
    r0 = norm(pos)
    s0 = dot(pos, vel)
    b  = 2gm/r0 - dot(vel, vel)

    # try to bracket the root
    # dt/dx = r >= 0, and monotonic, so we can step until we find a bracket
    # bracket = (0.0, dt/r)

    # better initial guesses
    xh = kepler_guess(pos, vel, dt, gm)
    # xh = dt/r0
    # xh = 1000.0
    yh, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh -= dt

    xl = zero(xh)
    # xl = 0.0
    yl = -dt
    rl = r0

    # establish the bracket
    i = 0
    while i < 1000 && (sign(yl) == sign(yh) || isinf(yh) || isnan(yh))
        i += 1
        if isinf(yh) || isnan(yh) #|| isnan(rh)
            xh = (xl + xh)/2
            yh, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh -= dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        elseif sign(yl) == sign(yh) && !isnan(rh)
            xl  = xh
            # xh += (dt - yh)/rh
            xh *= 2
            yl  = yh
            rl  = rh
            yh, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh -= dt
        end
    end

    return find_zero(_x -> universal_kepler(_x, b, r0, s0, gm) - dt, (xl, xh), A42())
end

test_func(x) = cos(x) - x

function kepler_guess(pos, vel, dt, gm)
    r0 = norm(pos)
    s0 = dot(pos, vel)
    b  = 2gm/r0 - dot(vel, vel)
    x0 = if abs(gm*b) < 1e-6
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
    else 
        throw((gm = gm, b = b))
        # dt/r0
    end
    return x0
end

function kepler_guess_canonical(pos, vel, dt)
    s0 = dot(pos, vel)
    b  = 2 - dot(vel, vel)
    x0 = if abs(b) < 1e-6
        # parabolic (Vallado)
        h = cross(pos, vel)
        p = dot(h, h)
        s = acot(3*sqrt(1/p^3)*dt)/2
        w = atan(cbrt(tan(s)))
        sqrt(p)*2*cot(2w)
    elseif b < 0
        # hyperbolic (Vallado)
        a = 1/b
        abs(sqrt(-a)*log(-2dt/(a*(s0+sqrt(-a)*(1 - 1/a)))))
    elseif b > 0
        # elliptic
        dt
    else 
        throw((gm = gm, b = b))
        # dt/r0
    end
    return x0
end

function solve_kepler_universal_lmm12(pos, vel, gm, dt; max_iter = 1000, tol = 1e-15)
    # constants
    r0 = norm(pos)
    s0 = dot(pos, vel)
    b  = 2gm/r0 - dot(vel, vel)

    # try to bracket the root
    # dt/dx = r >= 0, and monotonic, so we can step until we find a bracket
    # bracket = (0.0, dt/r)
    xl = 0.0
    yl = -dt
    rl = r0

    # better initial guesses
    xh = kepler_guess(pos, vel, dt, gm)
    # xh = dt/r0

    dth, rh = universal_kepler2(xh, b, r0, s0, gm)
    yh = dth - dt

    # establish the bracket
    i = 0
    while i < 1000 && (sign(yl) == sign(yh) || isinf(yh) || isnan(yh))
        i += 1
        if isinf(yh) || isnan(yh) #|| isnan(rh)
            xh = (xl + xh)/2
            # xh += (yh + dt)/rh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
            if xl == xh 
                throw("dt exceeds the computable range of values")
            end
        elseif sign(yl) == sign(yh) && !isnan(rh)
            xl  = xh
            # xh += (dth)/rh
            xh *= 2
            yl  = yh
            rl  = rh
            dth, rh = universal_kepler2(xh, b, r0, s0, gm)
            yh = dth - dt
        end
    end

    # we can now guaruntee a solution
    # TODO: utilize derivative based methods 
    # initialize at the end of the bracket, just for to have these defined for termination condition
    x = xh
    y = yh
    i = 0

    # safeguard against slow convergence by tracking iterations since either end of the bracket changed
    hi = 0
    li = 0
    n_force = 4
    # x1 = xl 
    # x2 = xh
    # while abs(xh - xl) > tol && i < max_iter && step > tol
    xb = (xl + xh)/2
    while abs(xh - xl) > tol && xb != xh && xb != xl && i < max_iter 
        i += 1
        if hi >= n_force || li >= n_force   
            # check how the brackets are progressing
            # if one hasnt changed in enough iterations, force a bisection step
            x = (xh + xl)/2
            
            # if a bisection returns a bracket end, we've converged
            if x == xh || x == xl
                break
            end
        else 
            x = lmm12_step(xl, xh, yl, yh, rl, rh)
        end
        
        y, r = universal_kepler2(x, b, r0, s0, gm)
        y   -= dt

        if y == 0
            break 
        end

        # assign the bracket
        if sign(y) == sign(yl)
            hi += 1
            li = 0
            xl = x
            yl = y
            rl = r
        elseif sign(y) == sign(yh)
            li += 1
            hi = 0
            xh = x
            yh = y
            rh = r
        else
            throw("bad step, sign check doesn't proeprly evaluate")
        end
    end

    if i == max_iter; throw("hit max iterations on kepler solve"); end

    return x
end

# function universal_kepler(y, l, k1, k2, k3)
#     _, c1, c2, c3 = stumpff(l*y^2)
#     L = y*(k1*c1 + y*(k2*c2 + y*k3*c3))
#     return L
# end

# function universal_kepler_param(x; b = 1.0, r0 = 1.0, s0 = 1.0, gm = 1.0, t = 1.0)
#     return universal_kepler(x, b, r0, s0, gm) - t
# end
# universal_kepler_param(x, b = 1.0, r0 = 1.0, s0 = 1.0, gm = 1.0, t = 1.0) = universal_kepler(x, b, r0, s0, gm) - t

function universal_kepler(x, b, r0, s0, gm)
    z  = b*x^2
    _, c1, c2, c3 = stumpff(z)
    # _, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    # dt = x*r0*c1 + (x^2)*s0*c2 + gm*(x^3)*c3
    return dt
end

function universal_kepler_canonical(x, b, s0)
    z  = b*x^2
    _, c1, c2, c3 = stumpff(z)
    dt = x*(c1 + x*(s0*c2 + x*c3))
    return dt
end

function universal_kepler2(x, b, r0, s0, gm)
    z  = b*x^2
    c0, c1, c2, c3 = stumpff(z)
    # c0, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    r  = r0*c0 + x*(s0*c1 + x*gm*c2)
    return dt, r
end


function universal_kepler2_canonical(x, b, s0)
    z  = b*x^2
    c0, c1, c2, c3 = stumpff(z)
    # c0, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(c1 + x*(s0*c2 + x*c3))
    r  = c0 + x*(s0*c1 + x*c2)
    return dt, r
end

function stumpff(z::T) where {T}
    if z > 1e-3
        sin2 = sin(sqrt(z)/2)
        cos2 = cos(sqrt(z)/2)

        c1 = 2sin2*cos2/sqrt(z)
        c2 = 2sin2^2/z
        c0 = 1 - z*c2
        c3 = (1 - c1)/z

        return c0, c1, c2, c3
    elseif z < -1e-3
        sin2 = sinh(sqrt(-z)/2)
        cos2 = cosh(sqrt(-z)/2)

        c1 = 2sin2*cos2/sqrt(-z)
        c2 = -2sin2^2/z
        c0 = 1 - z*c2
        c3 = (1 - c1)/z

        return c0, c1, c2, c3
    else
        # a1::T = 1/82
        # a2::T = 1/132
        # a3::T = 1/90
        # a4::T = 1/56
        # a5::T = 1/30
        # a6::T = 1/12
        # b1::T = 1/210
        # b2::T = 1/156
        # b3::T = 1/110
        # b4::T = 1/72
        # b5::T = 1/42
        # b6::T = 1/20
        # b7::T = 1/6
        # c2 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z*a1)*a2)*a3)*a4)*a5)*a6)/2
        # c3 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z*b1)*b2)*b3)*b4)*b5)*b6)*b7
        c3 = 1/6 - z/120 
        c2 = 1/2 - z/24
        p  = -z/720 # start with the 6! terms
        # i  = 6
        # while i < 15
        for i in 6:2:20
            p  *= z/i
            c2 += p
            p  /= i + 1
            c3 += p
        end
        c1 = 1 - z*c3
        c0 = 1 - z*c2
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
        c4 = 1/24  - z/720
        c5 = 1/120 - z/5040
        p  = -z/40320 # start with the 8! terms
        i  = 8
        # while i < 15
        for i in 8:2:20
            p  *= z/i
            c4 += p
            p  /= i + 1
            c5 += p
            i  += 2
        end
        c2 = 1/2 - z*c4
        c3 = 1/6 - z*c5
        c0 = 1 - z*c2
        c1 = 1 - z*c3
        return c0, c1, c2, c3, c4, c5
    end
end

function stumpff_fold(z)
    n  = 0
    zn = z
    while abs(zn) > 0.001
        zn /= 4
        n  += 1
    end
    # c2 = (1-zn*(1-zn*(1-zn*(1-zn*(1-zn*(1-zn*a1)*a2)*a3)*a4)*a5)*a6)/2
    # c3 = (1-zn*(1-zn*(1-zn*(1-zn*(1-zn*(1-zn*b1)*b2)*b3)*b4)*b5)*b6)*b7
    c5 = stumpff_continued(5, zn)
    c4 = stumpff_continued(4, zn)
    c3 = 1/6 - zn*c5
    c2 = 1/2 - zn*c4
    c1 = 1 - zn*c3
    while n > 0
        zn *= 4
        c5  = (c5 + c4 + c3*c2)/16
        c4  = c3*(1 + c1)/8
        c3  = 1/6 - zn*c5
        c2  = 1/2 - zn*c4
        c1  = 1 - zn*c3
        n  -= 1
    end
    c0 = 1 - zn*c2
    return c0, c1, c2, c3, c4, c5
end

""" evaluate the ith stumpff function as a continued fraction. Uses a fixed number of terms, so should be restricted to the regime near 0 (z < 1e-3 for double precision).
"""
function stumpff_continued(i, z)
    n = 6
    c = one(z)
    for k in n:-1:1
        # c -= z/((i + 2k)*(i + k))
        c = 1 - z/((i + 2k)*(i + k))
    end
    return c*factorial(i)
end

# """ 
# [DEPRECATED]
# Compute the position/velcoity along a Keplerian orbit after given timespan. Uses the universal variable formaulation
# from Danby, with initial guesses from Vallado, and modification to the stumpff functions from Wisdam + Hernandez 2015.

# normalizes the orbit to r0 = 1 and gm = 1, so there will likely be slight differences with other implementations from 
# floating point operations. This conveniently leads to s (Danby) equal to x (Vallado), and alpha (Vallado) equal to 
# beta (Wisdom, alpha in Danby). Total number of operations is reduced (though not likely impactful to most users)

# Uses a bracketing method to guaruntee convergence IF the upper bound guess is well behaved (the lower bound is 
# guarunteed). It is possible for the initial guess to evaluate to infinity if the timespan is too long, or if an orbit 
# is near enough to rectilinear that numerical issues lead to e = 1. While dropping something into the sun is perfectly
# well defined up to impact time, the singularity in the equations of motion causes a problem afterwards. More 
# importantly it invalidates the initial guess. There is an argument to be made that solve up to the point fof numeric 
# problems is desirable, (since we solve long time spans to the point of failure), and inward falling orbits are 
# physically reasonable, so future revisions may include such a solution. Until I find that solution, no rectilinear 
# orbits. 
# """
# function _propagate(pos, vel, dt, gm; max_iter = 20)
#     if all(vel .== 0) || all(pos .== 0) 
#         throw("invalid orbit")
#     end

#     if dt == 0
#         return (pos, vel)
#     elseif dt < 0 # THE STUPID GUESS DOESNT WORK FOR BACKWARDS TIME
#         posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
#         return posf, -velf
#     end

#     # normalize so that r0 = 1, gm = 1
#     # these are chosen to elliminate some computations
#     # and so that the numbers we're working with are (hopefully) in a nice floating point regime
#     DU   = norm(pos)
#     TU   = sqrt(DU^3 / gm)
#     pos0 = pos/DU
#     vel0 = vel/DU*TU
#     dt0  = dt/TU

#     dr0 = dot(vel0, pos0)
#     a   = 2.0 - dot(vel0, vel0) # 1/semimajor axis! alpha = beta in these coordinates

#     Hvec = cross(pos0, vel0)
#     Evec = cross(vel0, Hvec) - pos0
#     e    = norm(Evec)
    
#     # we have non-dimentionalized in a way that s = x, so we can mix equations from both Vallado and Danby
#     x = if abs(a) < 1e-6 
#         # near-parabolic
#         # converted from Vallado.
#         # need to study this
#         # NOT an upper bound, so still need to establish a bracket
#         p = dot(Hvec, Hvec)
#         s = 0.5acot(3dt*sqrt(1/p^3))
#         w = atan(tan(s)^(1/3))
#         2sqrt(p)*cot(2w)
#     else
#         x = dt0*a/(1 - e)
#     end

#     if x == Inf
#         throw("ERROR initial guess Inf. Likely causes are long timespans or near-rectilinear orbits. a = $(1/a*DU) e = $e q = $(a*(1-e)) dt = $dt pos = $pos vel = $vel gm = $gm")
#     end

#     # TODO: Add a check for near-rectilinear orbits. devise a solution

#     # set up the bracket
#     # only handle the forward time case (due to the above check)
#     x_br = 0.0
#     y_br = -dt0

#     # _, c1, c2, c3 = stumpff(a*x^2)
#     # y = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
#     _, g1, g2, g3 = stumpff(a, x)
#     y = g1 + dr0*g2 + g3 - dt0
    
#     while sign(y) == sign(y_br)
#         # take forward steps, maintaining the size of our interval
#         # we could use some newton steps here, and would probably work better
#         # but I'm lazy, and this is exceedingly likely to require at most 1 step
#         x_br = x
#         y_br = y
#         x    *= 2
#         # _, c1, c2, c3 = stumpff(a*x^2)
#         # y    = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
#         _, g1, g2, g3 = stumpff(a, x)
#         y = g1 + dr0*g2 + g3 - dt0
#     end

#     # we now have a bracket, and can find the root.
#     # chain a couple of closures together to get the right args to the right places
#     # TODO: implement a bracketed root finder using higher derivatives
#     bracket = (x_br, x)
#     # println("bracket = $bracket")
#     # _f1 = _x -> (_x, stumpff(a*_x^2)...)
#     # _f2 = ((_x, _, _c1, _c2, _c3),) -> _x*_c1 + dr0*_c2*_x^2 + _c3*_x^3
#     _f1 = _x -> stumpff(a, _x)
#     _f2 = ((_, _g1, _g2, _g3),) -> _g1 + dr0*_g2 + _g3
#     x = find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, A42())

#     # _, c1, c2, c3 = stumpff(a*x^2)
#     g0, g1, g2, _ = stumpff(a, x)
#     f    = 1.0 - g2
#     g    = g1 + dr0*g2
#     posf = f*pos0 + g*vel0

#     rf = norm(posf)
#     df = -g1/rf
#     dg = (g0 + dr0*g1)/rf
#     # dg   = 1.0 - g2/rf
#     # df   = x*(a*c3*x^2 - 1.0)/rf
#     velf = df*pos0 + dg*vel0

#     # re-dimensionalize our inputs
#     return posf*DU, velf*DU/TU
# end