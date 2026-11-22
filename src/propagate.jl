""" partitioned Cartesian State Transition Matrix. 
"""
struct KeplerSTM{M}
    dX_dX0::M
    dX_dV0::M
    dV_dX0::M
    dV_dV0::M
end

""" Takes advantage of the symplectic structure of the STM to quickly invert 

Note that the inverted STM is also the time reversed STM
"""
function Base.inv(stm::KeplerSTM)
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

function _propagate_with_partials(state::Cartesian, t)
    posf, velf, stm = propagate_with_partials(state.position, state.velocity, t - state.epoch, state.gm)
    return Cartesian(posf, velf, t, state.gm), stm
end

# NOTE: THE f used here is 1 - f', where f' is the traditional f
function propagate(pos, vel, dt, gm)
    delp, delv = kepler_step(pos, vel, dt, gm)
    return pos + delp, vel + delv
end

function propagate_with_partials(pos, vel, dt, gm)
    delp, delv, stm = kepler_step_with_stm(pos, vel, dt, gm)
    return pos + delp, vel + delv, stm
end

function propagate_old(pos, vel, dt, gm)
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
function propagate_with_partials_old(pos, vel, dt, gm)
    DU = norm(pos)*sign(gm)
    TU = sqrt(DU^3/gm)

    pos /= DU
    vel /= DU/TU
    dt  /= TU

    # constants
    b  = 2 - dot(vel, vel)
    s0 = dot(vel, pos)

    # x = solve_kepler_universal_normalized_new(pos, vel, dt)
    x = solve_universal_kepler(dt, 1.0, s0, b, 1.0)

    # compute f and g functions
    # _, c1, c2, c3, c4, c5 = stumpff5(b*x^2)
    U0, U1, U2, _, U4, U5 = universal05(b, x)
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

    # compute the partials (Goodyear)
    U  = U2*dt + (x*U4 - 3*U5)
    X0 = transpose(hcat(pos, vel))
    Xf = hcat(posf, velf)
    acc0 = -normalize(pos)
    accf = -posf/rf^3
    
    M1 = transpose(SMatrix{2, 2}(
        -(df*U1 + f),
        -df*U2,
        f*U1,
        f*U2,
    ))
    M2 = transpose(SMatrix{2, 2}(-df, -dg, f, g))*U2
    M3 = transpose(SMatrix{2, 2}(
        -df*(U0 + 1/rf + rf)/rf,
        -(df*U1 + dg/rf)/rf,
        df*U1 + f,
        df*U2
    ))
    M4 = transpose(SMatrix{2, 2}(
        -(df*U1 + dg/rf)/rf,
        -dg*U1/rf,
        df*U2,
        dg*U2,
    ))

    dxdx =  (
            (1.0 + f)*I
        + U*velf*transpose(acc0) 
        # + velf*transpose(acc0)*DU^2/TU^3
        + Xf*M1*X0
    )
    dxdv =          g*I - U*velf*transpose(vel)  + Xf*M2*X0
    dvdx =         df*I + U*accf*transpose(acc0) + Xf*M3*X0
    dvdv = (1.0 + dg)*I - U*accf*transpose(vel)  + Xf*M4*X0
    
    # return posf*DU, velf*DU/TU, dxdx, dxdv*TU, dvdx/TU, dvdv
    return posf*DU, velf*DU/TU, KeplerSTM(dxdx, dxdv*TU, dvdx/TU, dvdv)
end

# TODO: Add non-normalized variants
# NOTE: THE f used here is f' - 1, where f' is the traditional f (as used in Battin)
#     : likewise for dg 
# NOTE: normalized to r0 = 1, gm = 1
# stm_pos_pos0_normalized(pos, posf, vel, velf, rf, f, C) = (
#     (1.0 + f)*I
#     + rf*(velf - vel)*transpose((velf - vel)) 
#     - f*posf*transpose(pos) 
#     + C*velf*transpose(pos)
# )
# 
# stm_pos_vel0_normalized(pos, posf, vel, velf, f, g, C) = (
#     g*I 
#     - f*(
#           (posf - pos)*transpose(vel) 
#         - (velf - vel)*transpose(pos)
#         )
#     + C*velf*transpose(vel)
# )
# 
# stm_vel_pos0_normalized(pos, posf, vel, velf, rf, df, C) = (
#     - (velf - vel)*transpose(pos)
#     - (1/rf^2)*posf*transpose(velf - vel) 
#     + df*(I - (1/rf^2)*posf*transpose(posf) + (1/rf)*(posf*transpose(velf) - velf*transpose(posf))*posf*transpose(velf - vel))
#     - C/(rf)^3*posf*transpose(pos)
# )
# 
# stm_vel_vel0_normalized(pos, posf, vel, velf, rf, f, dg, C) = (
#       (1.0 + dg)*I
#     + (velf - vel)*transpose((velf - vel)) 
#     + (1/rf^3)*(
#         - f*posf*transpose(pos)
#         - C*posf*transpose(vel)
#     )
# )

function kepler_step(pos, vel, dt, gm)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)

    s = solve_universal_kepler(dt, r0, rdr0, b, gm)
    f, g, df, dg = lagrange_step_coefficients(pos, vel, gm, s)
    delp =  f*pos +  g*vel
    delv = df*pos + dg*vel
    return delp, delv
end

function kepler_step_with_stm(pos, vel, dt, gm)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)

    s = solve_universal_kepler(dt, r0, rdr0, b, gm)
    f, g, df, dg, stm = lagrange_step_coefficients_with_stm(pos, vel, dt, gm, s)
    delp =  f*pos +  g*vel
    delv = df*pos + dg*vel
    return delp, delv, stm
end

function lagrange_step_coefficients(pos, vel, gm, s)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)
    
    U = universal03(b, s) # index i -> i + 1
    f = -gm/r0*U[3]
    g = r0*U[2] + rdr0*U[3]
    
    posf = pos + f*pos + g*vel
    rf   = norm(posf)

    df = -gm/(r0*rf)*U[2]
    dg = -gm/rf*U[3]
    return f, g, df, dg
end

function lagrange_step_coefficients_with_stm(pos, vel, dt, gm, s)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)
    
    U0, U1, U2, _, U4, U5 = universal05(b, s)
    f, g, df, dg = lagrange_step_coefficients(pos, vel, gm, s)

    posf = pos + f*pos + g*vel
    velf = vel + df*pos + dg*vel

    rf = norm(posf)

    acc0 = -gm*pos  / r0^3
    accf = -gm*posf / rf^3

    x0 = transpose(hcat(pos, vel))
    xf = hcat(posf, velf)

    # from Goodyear 1965 + errata
    C = dt*U2 + gm*(s*U4 - 3*U5) # U from Goodyear
    M1 = SMatrix{2, 2}(
        -(df*U1 + f/r0)/r0,
        -df*U2,
        f*U1/r0,
        f*U2
    )
    M2 = SMatrix{2, 2}(-df, -dg, f, g)*U2
    M3 = SMatrix{2, 2}(
        -df*(U0 + r0/rf + rf/r0)/r0/rf,
        -(df*U1 + dg/rf)/rf,
        (df*U1 + f/r0)/r0,
        df*U2
    )
    M4 = SMatrix{2, 2}(
        -(df*U1 + dg/rf)/rf,
        -dg*U1/rf,
        df*U2,
        dg*U2
    )

    STM = KeplerSTM(
        (
            I*(1 + f)
            + C*velf*transpose(acc0)
            # + velf*transpose(acc0)
            + xf*transpose(M1)*x0
        ),
        I*g        - C*velf*transpose(vel)  + xf*transpose(M2)*x0,
        I*df       + C*accf*transpose(acc0) + xf*transpose(M3)*x0,
        I*(1 + dg) - C*accf*transpose(vel)  + xf*transpose(M4)*x0,
    )
    
    return f, g, df, dg, STM
end

function solve_universal_kepler(dt, r0, rdr0, b, gm)
    # establish initial guess using a newton step from zero (what about a halley step?)
    sl = zero(dt)
    sh = 1.01*dt/r0 # deliberately try to over shoot too good of a guess

    # make sure it's a bracket
    fl = -dt
    fh = universal_kepler(r0, rdr0, b, gm, sh) - dt
    while sign(fl)*sign(fh) > 0 || (isfinite(fl) && !isfinite(fh))
        if isfinite(fh)
            # undershot, shift the bracket
            (sl, fl) = (sh, fh)
            sh *= 2
        else
            # overflow but lower bound still valid, bisect
            sh = sl + (sh - sl)/2
        end
        fh = universal_kepler(r0, rdr0, b, gm, sh) - dt
    end

    if fl == 0
        return sl
    end

    if fh == 0
        return sh
    end

    @assert isfinite(fl) "Failed to find a finite bracket given input type $(typeof(fh))"
    
    # return find_zero(s -> universal_kepler(r0, rdr0, b, gm, s) - dt, (sl, sh), Roots.Chandrupatla())
    return chandrupatla(s -> universal_kepler(r0, rdr0, b, gm, s) - dt, (sl, sh))
    # function univ_kep2(s)
    #     t, r = universal_kepler_with_deriv(r0, rdr0, b, gm, s)
    #     return (t - dt, r)
    # end
    # return newton_bracket(univ_kep2, (sl, sh))
end

function universal_kepler(r0, rdr0, b, gm, s)
    U = universal03(b, s)
    t = r0*U[2] + rdr0*U[3] + gm*U[4]
    return t
end

function universal_kepler_with_deriv(r0, rdr0, b, gm, s)
    U = universal03(b, s)
    t = r0*U[2] + rdr0*U[3] + gm*U[4]
    r = r0*U[1] + rdr0*U[2] + gm*U[3]
    return t, r
end
