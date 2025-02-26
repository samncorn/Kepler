module Kepler

using LinearAlgebra
using Roots
using Printf

# Write your package code here.
include("root-finding.jl")
include("stumpff.jl")

""" 
Keplerian orbit evolution. Uses the universal variable formaulation from Danby, with initial guesses from Vallado, and
modification to the stumpff functions from Wisdam + Hernandez 2015
"""
function solve(pos0, vel0, dt, gm; max_iter = 20, anomaly_tol = 1e-8, parabolic_tol = 1e-12)
    if dt == 0
        return (pos0, vel0)
    end

    r0  = norm(pos0)
    dr0 = dot(vel0, pos0) / r0
    v02 = dot(vel0, vel0)
    alpha = 2gm/r0 - v02

    # solve the cubic
    # what if multiple real roots?
    s = begin # if abs(alpha) > parabolic_tol # elliptic or hyperbolic, use the simple bracketing 
        Hvec = cross(pos0, vel0)
        Evec = cross(vel0, Hvec)/gm - pos0/r0
        e = norm(Evec)
        dt*alpha/(gm*(1 - e))
    end
    # else # 
    #     a0 = -6dt/gm
    #     a1 = 6r0/gm
    #     a2 = 3r0*dr0/gm

    #     q = a1/3 - (a2^2)/9
    #     r = (a1*a2 - 3a0)/6 - (a2^3)/27
    #     p1 = (r + sqrt(q^3 + r^2))^(1/3)
    #     p2 = (r - sqrt(q^3 + r^2))^(1/3)
    #     p1 + p2 - a2/3
    # end
    # set up the bracket
    br = 0.0
    br_y = universal_kepler(br, alpha, r0, dr0, gm, dt)
    y    = universal_kepler(s, alpha, r0, dr0, gm, dt)
    while sign(y) == sign(br_y)
        br = s
        br_y = y
        s = 2s
        y = universal_kepler(s, alpha, r0, dr0, gm, dt)
    end
    bracket = (min(br, s), max(br, s))
    s = find_zero(s -> universal_kepler(s, alpha, r0, dr0, gm, dt), bracket, A42())

    # n_iter = 1
    # while !converged
    #     y, dy, ddy = universal_kepler3(s, alpha, r0, dr0, k, dt)
    #     ds1 = y / dy
    #     ds2 = y / (dy - (ds1 * ddy) / 2)

    #     ds = ds2
    #     if sign(-ds) != sign(s - br) # bracket increasing, take a secant step

    #     end

    #     # new bracket

    #     converged = abs(s - br) < anomaly_tol || n_iter == max_iter
    #     n_iter += 1
    # end


    # f1 = x -> universal_kepler3(x, alpha, r0, dr0, gm)
    # f2 = x -> (x[1] - dt, x[2:end]...)
    # f  = x -> f2(f1(x))
    # s, n_iter = laguerre(f, s0, anomaly_tol, max_iter)

    # if n_iter >= max_iter
    #     H = cross(pos0, vel0)
    #     E = cross(vel0, H)/gm - pos0/r0
    #     ecc = norm(E)
    #     @warn @sprintf "kepler solverer failed to converge after %i iterations" n_iter
    #     @warn @sprintf "position: %.7e %.7e %.7e" pos0...
    #     @warn @sprintf "velocity: %.7e %.7e %.7e" vel0...
    #     @warn @sprintf "gm: %.7e" gm
    #     @warn @sprintf "dt: %.7e" dt
    #     @warn @sprintf "α : %.7e" alpha
    #     @warn @sprintf "e : %.7e" ecc
    #     return NaN*pos0, NaN*vel0
    # end

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

function universal_kepler(s, alpha, r0, dr0, k, dt)
    _, c1, c2, c3 = stumpff(alpha*s^2)
    d = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3 - dt
    return d
end

function universal_kepler2(s, alpha, r0, dr0, k, dt)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    d = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3 - dt
    d1 = r0*c0 + r0*dr0*s*c1 + k*s^2*c2
    return d, d1
end

function universal_kepler3(s, alpha, r0, dr0, k, dt)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    d = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3 - dt
    d1 = r0*c0 + r0*dr0*s*c1 + k*s^2*c2
    # d2 = (-r0*alpha + k)*s*c1 + r0*dr0*c0
    d2 = r0*dr0*c0 + k*(1 - alpha*r0)*s*c1
    return d, d1, d2
end

function universal_kepler4(s, alpha, r0, dr0, k, dt)
    c0, c1, c2, c3 = stumpff(alpha*s^2)
    d = r0*s*c1 + r0*dr0*s^2*c2 + k*s^3*c3 - dt
    d1 = r0*c0 + r0*dr0*s*c1 + k*s^2*c2
    # d2 = (-r0*alpha + k)*s*c1 + r0*dr0*c0
    d2 = r0*dr0*c0 + k*(1 - alpha*r0)*s*c1
    # d3 = (-r0*alpha + k)*c0 - r0*dr0*alpha*s*c1
    d3 = k*((1 - alpha*r0)*c0 - alpha*r0*dr0*s*c1)
    return d, d1, d2, d3
end

end # module