""" partitioned Cartesian State Transition Matrix. 
"""
struct STM{M}
    dX_dX0::M
    dX_dV0::M
    dV_dX0::M
    dV_dV0::M
end

""" Takes advantage of the symplectic structure of the STM to quickly invert 

Note that the inverted STM is also the time reversed STM
"""
function Base.inv(stm::STM)
    return STM(
         transpose(stm.dV_dV0),
        -transpose(stm.dX_dV0),
        -transpose(stm.dV_dX0),
         transpose(stm.dX_dX0),
    )
end

"Universal kepler solver."
function propagate(state::Cartesian, t)
    posf, velf = propagate(state.position, state.velocity, t - state.epoch, state.gm)
    return Cartesian(posf, velf, t, state.gm)
end

function propagate_stm(state::Cartesian, t)
    posf, velf, stm = propagate_stm(state.position, state.velocity, t - state.epoch, state.gm)
    return Cartesian(posf, velf, t, state.gm), stm
end

# adopt a autodiff like interface?
# propagate(pos, )

# NOTE: THE f used here is 1 - f', where f' is the traditional f
function propagate(pos, vel, dt, gm)
    DU = norm(pos)
    TU = sqrt(DU^3/abs(gm))
    pos /= DU
    vel /= DU/TU
    dt  /= TU

    b  = 2.0 - dot(vel, vel)
    s0 = dot(vel, pos)

    x = solve_kepler_universal_normalized_new(pos, vel, dt)

    # compute f and g functions
    _, U1, U2, _ = universal03(b, x)

    # we use the modified versions of f and dg given by Rein et al
    f = -U2
    g = U1 + s0*U2
    posf = f*pos + g*vel + pos
    
    rf = norm(posf)
    df = -U1/rf
    dg = -U2/rf
    velf = df*pos + dg*vel + vel

    return posf*DU, velf*DU/TU
end

""" Returns both the propagated state and the Keplerian STM.

See Battin 9.7
"""
function propagate_stm(pos, vel, dt, gm)
    if dt == 0
        return pos, vel
    end

    DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    # constants
    b  = 2 - dot(vel, vel)
    s0 = dot(vel, pos)

    x = solve_kepler_universal_normalized_new(pos, vel, dt)

    # compute f and g functions
    # _, c1, c2, c3, c4, c5 = stumpff5(b*x^2)
    _, U1, U2, _, U4, U5 = universal05(b, x)
    # U1 = x*c1
    # U2 = x^2*c2
    # U4 = x^4*c4
    # U5 = x^5*c5

    # we use the modified versions of f and dg given by Rein et al
    f = -U2
    g = U1 + s0*U2
    posf = f*pos + g*vel + pos
    
    rf = norm(posf)
    df = -U1/rf
    dg = -U2/rf
    velf = df*pos + dg*vel + vel

    # compute the partials (Battin)
    # C  = (x^2)*((x^3)*(3c5 - c4) - dt*c2)
    C    = 3U5 - x*U4 - dt*U2
    dxdx = stm_pos_pos0_normalized(pos, posf, vel, velf, rf, f, C)     # units of DU/DU = 1
    dxdv = stm_pos_vel0_normalized(pos, posf, vel, velf, f, g, C)*TU   # Units of TU
    dvdx = stm_vel_pos0_normalized(pos, posf, vel, velf, rf, df, C)/TU # Units of 1/TU
    dvdv = stm_vel_vel0_normalized(pos, posf, vel, velf, rf, f, dg, C) # Units of DU/TU*TU/DU = 1

    # return posf*DU, velf*DU/TU, dxdx, dxdv*TU, dvdx/TU, dvdv
    return posf*DU, velf*DU/TU, STM(dxdx, dxdv, dvdx, dvdv)
end

# TODO: Add non-normalized variants
# NOTE: THE f used here is f' - 1, where f' is the traditional f (as used in Battin)
#     : likewise for dg 
# NOTE: normalized to r0 = 1, gm = 1
stm_pos_pos0_normalized(pos, posf, vel, velf, rf, f, C) = (
    (1.0 + f)*I
    + rf*(velf - vel)*transpose((velf - vel)) 
    - f*posf*transpose(pos) 
    + C*velf*transpose(pos)
)

stm_pos_vel0_normalized(pos, posf, vel, velf, f, g, C) = (
    g*I 
    - f*(
          (posf - pos)*transpose(vel) 
        - (velf - vel)*transpose(pos)
        )
    + C*velf*transpose(vel)
)

stm_vel_pos0_normalized(pos, posf, vel, velf, rf, df, C) = (
    - (velf - vel)*transpose(pos)
    - (1/rf^2)*posf*transpose(velf - vel) 
    + df*(I - (1/rf^2)*posf*transpose(posf) + (1/rf)*(posf*transpose(velf) - velf*transpose(posf))*posf*transpose(velf - vel))
    - C/(rf)^3*posf*transpose(pos)
)

stm_vel_vel0_normalized(pos, posf, vel, velf, rf, f, dg, C) = (
      (1.0 + dg)*I
    + (velf - vel)*transpose((velf - vel)) 
    + (1/rf^3)*(
        - f*posf*transpose(pos)
        - C*posf*transpose(vel)
    )
)

function solve_kepler_universal_normalized_new(pos, vel, dt)
    s0 = dot(vel, pos)
    b  = 2.0 - dot(vel, vel)

    # use x = 0 as the initial contrapoint
    x1 = zero(dt)
    y1 = -dt
    r1 = one(dt)
    # s1 = s0
    p1 = (x = x1, y = y1, dy = r1)

    # x2     = kepler_guess_canonical(pos, vel, dt)
    x2     = dt
    # x2 = -dt/(1.0 - dt*s0/2)
    # x2 = dt + (s0/2)*dt^2
    y2, r2 = universal_kepler2_canonical(x2, b, s0)
    y2    -= dt
    p2 = (x = x2, y = y2, dy = r2)

    # bracket
    while sign(p2.y) == sign(p1.y) && isfinite(p1.y)
        if p2.y == Inf
            # println("halving")
            x2 /= 2
            y2, r2 = universal_kepler2_canonical(x2, b, s0)
            y2    -= dt
            p2 = (x = x2, y = y2, dy = r2)
        else
            # println("stepping")
            # try to bracket
            # compute a newton step, and deliberately overshoot (double the step)
            p1 = p2

            # x2 = (x2*r2 - 2y2)/r2
            x2 *= 2
            y2, r2 = universal_kepler2_canonical(x2, b, s0)
            y2    -= dt
            p2 = (x = x2, y = y2, dy = r2)
        end
    end

    if y2 == 0
        return x2
    end

    if y2 == Inf || isnan(y2)
        throw("OVERFLOW")
    end

    # one manual inverse interpolation step
    x = flmsm1_step(p1, p2)
    # check for leaving the interval 
    if x < p1.x || x > p2.x
        x = (p1.x + p2.x)/2
    end
    y, r = universal_kepler2_canonical(x, b, s0)
    y   -= dt
    p3   = (x = x, y = y, dy = r)

    if y == 0
        return x
    end

    # shift the bracket, so that p1 is always the contrapoint
    # that is, p2 and p3 bracket the root
    if sign(p3.y) == sign(p2.y)
        (p1, p2) = (p2, p1)
    end

    # inverse interpolation until convergence
    # TODO: Try to prove convergence?
    # TODO: 2nd derivative methods
    i = 0
    # x = NaN # just to initialize the loop
    while p1.x != x && p2.x != x && i < 1_000
    # while abs(y) > 1e-12 && i < 1_000
        i += 1

        signab = sign((p1.y - p2.y)/(p1.x - p2.x))
        if signab == sign(p1.dy) && signab == sign(p2.dy) && signab == sign(p3.dy) && p1.y != p3.y
            # derivatives are well suited for inverse interpolation
            # if p2.y != p3.y
            #     # need unique values
            #     x = flmsm1_step(p1, p2, p3)
            # else
            #     x = flmsm1_step(p1, p3)
            # end

            # there is a numeric instability in the 3-point method that leads to a biased drift in some classes of orbits
            # the 2-point method seems to not encounter this, with only a modest hit to the convergence rate (2.73 vs 2.91)
            # so it's likely not essential to 
            # x = flmsm1_step(p1, p3)
            x = flmsm1_step(p1, p3)
        else
            # attempt newton's method
            x = (p3.x*p3.dy - p3.y)/p3.dy
        end

        # if we have stepped outside the bracket, bisect
        if (x < p3.x && x < p2.x) || (x > p3.x && x > p2.x)
            # x = (p2.x + p3.x)/2
            x = p2.x + (p3.x - p2.x)/2
        end

        if isnan(x)
            throw("NAN x1 = $(p1.x) x2 = $(p2.x) y1 = $(p1.y) y2 = $(p2.y) x3 = $(p3.x) y3 = $(p3.y)")
        end
        y, r = universal_kepler2_canonical(x, b, s0)
        y   -= dt
        p    = (x = x, y = y, dy = r)

        if y == 0
            return x
        end

        p1 = p2
        p2 = p3
        p3 = p
        if sign(p3.y) == sign(p2.y)
            (p1, p2) = (p2, p1)
        end
    end

    if i == 1000
        throw("Kepler solve hit max iterations")
    end

    return x
end

function solve_kepler_universal_normalized(pos, vel, dt)
    # r0 = dot(pos, pos)
    # @assert abs(1 - r0) < 

    # constants
    s0 = dot(vel, pos)
    b  = 2 - dot(vel, vel)

    # better initial guesses?
    # xh = kepler_guess_canonical(pos, vel, dt)
    xh = dt
    yh, rh = universal_kepler2_canonical(xh, b, s0)
    yh -= dt

    if yh == 0
        return xh
    end

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
        else 
            throw("WARNING SOMETHING WENT HORRIBLY WRONG")
        end
    end

    # return find_zero(_x -> universal_kepler_canonical(_x, b, s0) - dt, (xl, xh), A42())
    return chandrupatla_brent(_x -> universal_kepler_canonical(_x, b, s0) - dt, sort((xl, xh)))
    # return inverse_quintic(_x -> universal_kepler2_canonical_offset(_x, b, s0, dt), (xl, xh))
end

function solve_kepler_universal(pos, vel, gm, dt)
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
        else
            throw("WARNING SOMETHING WENT HORRIBLY WRONG")
        end
    end

    return find_zero(_x -> universal_kepler(_x, b, r0, s0, gm) - dt, (xl, xh), A42())
end

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
    b  = 2.0 - dot(vel, vel)
    x0 = if b > 1e-6
        # elliptic
        dt
    elseif b < -1e-6
        # hyperbolic
        a  = 1/b
        ca = sqrt(-a)
        sign(dt)*ca*log(-2dt*b/(s0 + sign(dt)*ca*(1.0 - b)))
    elseif abs(b) < 1e-6
        # parabolic
        h = cross(pos, vel)
        p = dot(h, h)
        s = acot(3*sqrt(1/p^3)*dt)/2
        w = atan(cbrt(tan(s)))
        sqrt(p)*2*cot(2w)
    else
        throw("invalid input b = $b")
    end
    return x0
end

function universal_kepler(x, b, r0, s0, gm)
    z  = b*x^2
    _, c1, c2, c3 = stumpff(z)
    # _, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    # dt = x*r0*c1 + (x^2)*s0*c2 + gm*(x^3)*c3
    return dt
end

# function universal_kepler_canonical(x, b, s0)
#     z  = b*x^2
#     _, c1, c2, c3 = stumpff(z)
#     dt = x*(c1 + x*(s0*c2 + x*c3))
#     return dt
# end

function universal_kepler2(x, b, r0, s0, gm)
    z  = b*x^2
    c0, c1, c2, c3 = stumpff(z)
    # c0, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(r0*c1 + x*(s0*c2 + gm*x*c3))
    r  = r0*c0 + x*(s0*c1 + x*gm*c2)
    return dt, r
end

# function universal_kepler2_canonical(x, b, s0)
#     z  = b*x^2
#     c0, c1, c2, c3 = stumpff(z)
#     # c0, c1, c2, c3, _, _ = stumpff_fold(z)
#     dt = x*(c1 + x*(s0*c2 + x*c3))
#     r = c0 + x*(s0*c1 + x*c2)
#     return dt, r
# end

function universal_kepler_canonical(x, b, s0)
    _, U1, U2, U3 = universal03(b, x)
    dt = U1 + s0*U2 + U3
    return dt
end

function universal_kepler2_canonical(x, b, s0)
    U0, U1, U2, U3 = universal03(b, x)
    dt = U1 + s0*U2 + U3
    r  = U0 + s0*U1 + U2 
    return dt, r
end

function universal_kepler2_canonical_offset(x, b, s0, _dt)
    z  = b*x^2
    c0, c1, c2, c3 = stumpff(z)
    # c0, c1, c2, c3, _, _ = stumpff_fold(z)
    dt = x*(c1 + x*(s0*c2 + x*c3))
    r = c0 + x*(s0*c1 + x*c2)
    return dt - _dt, r
end