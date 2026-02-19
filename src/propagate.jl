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
    posf, velf = propagate(state.position, state.velcoity, t - state.epoch, state.gm)
    return Cartesian(posf, velf, t, state.gm)
end

function propagate_stm(state::Cartesian, t)
    posf, velf, stm = propagate_stm(state.position, state.velcoity, t - state.epoch, state.gm)
    return Cartesian(posf, velf, t, state.gm), stm
end

# adopt a autodiff like interface?
# propagate(pos, )

# NOTE: THE f used here is 1 - f', where f' is the traditional f
function propagate(pos, vel, dt, gm)
    if dt == 0
        return pos, vel
    end

    if dt < 0
        posf, velf = propagate(pos, -vel, -dt, gm)
        return posf, -velf
    end

    DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    b = 2.0 - dot(vel, vel)

    x = solve_kepler_universal_A42_canonical(pos, vel, dt)

    # compute f and g functions
    _, c1, c2, c3 = stumpff(b*x^2)

    # more numerically precise? According to Rein + Tamayo WHFast paper
    f = -c2*x^2
    g = dt - c3*x^3
    posf = f*pos + g*vel + pos
    
    rf = norm(posf)
    df = -x*c1/rf
    dg = -(c2/rf)*x^2
    velf = df*pos + dg*vel + vel

    return posf*DU, velf*DU/TU
end

""" Returns both the propagated state and the Keplerian STM.

See Battin 9.7
"""
function propagate_stm(pos, vel, dt, gm)
    if dt == 0
        return pos, vel, I3, I3, I3, I3
    end

    if dt < 0
        # posf, velf = propagate(pos, -vel, -dt, gm; max_iter = max_iter)
        posf, velf, dxdx, dxdv, dvdx, dvdv = propagate_with_partials(pos, -vel, -dt, gm)
        return posf, -velf, dxdx, -dxdv, -dvdx, dvdv
    end

    DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    # constants
    b = 2 - dot(vel, vel)

    x = solve_kepler_universal_A42_canonical(pos, vel, dt)

    # compute f and g functions
    _, c1, c2, c3, c4, c5 = stumpff5(b*x^2)

    # more numerically precise? According to Rein + Tamayo WHFast paper
    f = -c2*x^2
    g = dt - c3*x^3
    posf = f*pos + g*vel + pos
    
    rf = norm(posf)
    df = -x*c1/rf
    dg = -(c2/rf)*x^2
    velf = df*pos + dg*vel + vel

    # compute the partials (Battin)
    C  = (x^2)*((x^3)*(3c5 - c4) - dt*c2)
    dxdx = stm_pos_pos0_normalized(pos, posf, vel, velf, rf, f, C)     # units of DU/DU = 1
    dxdv = stm_pos_vel0_normalized(pos, posf, vel, velf, f, g, C)*TU   # Units of TU
    dvdx = stm_vel_pos0_normalized(pos, posf, vel, velf, rf, df, C)/TU # Units of 1/TU
    dvdv = stm_vel_vel0_normalized(pos, posf, vel, velf, rf, f, dg, C) # Units of DU/TU*TU/DU = 1

    # return posf*DU, velf*DU/TU, dxdx, dxdv*TU, dvdx/TU, dvdv
    return posf*DU, velf*DU/TU, STM(dxdx, dxdv, dvdx, dvdv)
end

# TODO: Add non-normalized variants
# NOTE: THE f used here is 1 - f', where f' is the traditional f (as used in Battin)
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
        - (vel - velf)*transpose(pos)
        )
    + C*velf*transpose(vel)
)

stm_vel_pos0_normalized(pos, posf, vel, velf, rf, df, C) = (
    - (velf - vel)*transpose(pos)
    - (1/rf^2)*posf*transpose(velf - vel) 
    + df*(
        I - (1/rf^2)*posf*transpose(posf) 
          + (1/rf)*posf*transpose(delv)*(
              posf*transpose(velf)
            - velf*transpose(posf)
          )
        )
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

function solve_kepler_universal_canonical(pos, vel, dt)
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
        else 
            throw("WARNING SOMETHING WENT HORRIBLY WRONG")
        end
    end

    return find_zero(_x -> universal_kepler_canonical(_x, b, s0) - dt, (xl, xh), A42())
    # return chandrupatla_brent(_x -> universal_kepler_canonical(_x, b, s0) - dt, (xl, xh))
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
    r = c0 + x*(s0*c1 + x*c2)
    return dt, r
end
