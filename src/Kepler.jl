module Kepler

using LinearAlgebra

# Write your package code here.
include("root-finding.jl")
include("stumpff.jl")

""" 
"""
function propagate(pos0, vel0, dt, gm; max_iter = 20, anomaly_tol = 1e-8, parabolic_tol = 1e-6)
    if dt == 0
        return (pos0, vel0)
    end

    r0  = norm(pos0)
    dr0 = dot(vel0, normalize(pos0))
    v02 = dot(vel0, vel0)
    alpha = 2/r0 - v02/gm

    # initial guesses
    s0 = if alpha > parabolic_tol # should tune this param
        # elliptical
        dt*alpha*sqrt(gm)
    elseif alpha < -parabolic_tol
        # hyperbolic
        a = 1/alpha
        sign(dt)*sqrt(-a)*log(-2*gm*alpha*dt/(r0*dr0 + sign(dt)*sqrt(-gm*a)*(1 - r0*alpha)))
    else
        h = cross(pos0, vel0)
        p = dot(h, h) / gm
        s = 0.5acot(3sqrt(gm/p^3)*dt)
        w = atan((tan(s))^(1/3))
        2sqrt(p)*cot(2w)
    end

    f1 = x -> universal_kepler2(x, alpha, r0, dr0, gm)
    # f1 = x -> universal_kepler3(b, x, r0, dr0, gm)
    f2 = x -> (x[1] - dt, x[2:end]...)
    f = x -> f2(f1(x))
    s, n_iter = newton2(f, s0, anomaly_tol, max_iter)

    if n_iter >= max_iter
        @warn "kepler solverer failed to converge"
        return NaN*pos0, NaN*vel0
    end

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

function universal_kepler(s, alpha, r0, dr0, k)
    _, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3
    return dt
end

function universal_kepler2(s, alpha, r0, dr0, k)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3
    d1 = r0*c0 + r0*dr0*s*c1 + k*s^2*c2
    return dt, d1
end

function universal_kepler3(s, alpha, r0, dr0, k)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3
    d1 = r0*c0 + r0*dr0*s*c1 + k*s^2*c2
    d2 = (-r0*alpha + gm)*s*c1 + r0*dr0*c0
    return dt, d1, d2
end

function universal_kepler4(s, alpha, r0, dr0, k)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    dt = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3
    d1 = r0*c0 + r0*dr0*c1 + k*c2
    d2 = (-r0*alpha + gm)*s*c1 + r0*dr0*c0
    d3 = (-r0*alpha + gm)*c0 - r0*dr0*alpha*s*c1
    return dt, d1, d2, d3
end

end # module