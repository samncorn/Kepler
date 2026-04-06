""" Perform IOD using Herget's method. Assumes less than a single revolution.

The reference observations may be duplicates of the entries in list observations, to aid in sampling different reference pairs.

In the case of only 3 measurements, this method is very similar to a method laid out by Gooding, though with a slightly different formulation
"""
function herget_iod(observations, obs1, obs2, rho1, rho2, gm, c)
    x0 = SVector{2}(rho1, rho2)
    (rho1f, rho2f), _, _ = Kepler.least_squares(x -> herget_kernel(x, observations, obs1, obs2, gm, c), x0)

    return herget_solve(obs1, obs2, rho1f, rho2f, gm, c)
end

""" returns an iterator over the observations, each returning a tuple of ((dra, ddec), H)
"""
function herget_kernel(rho12, observations, obs1, obs2, gm, c; del = 0.001)
    return Iterators.map(o -> herget_residuals_with_partials(o, obs1, obs2, rho12[1], rho12[2], gm, c; del = del), observations)
end

# assumes 0-rev lambert solution
function herget_solve(obs1, obs2, rho1, rho2, gm, c)
    pos1 = rho1*obs1.angles + obs1.position
    pos2 = rho2*obs2.angles + obs2.position

    t1 = obs1.time - rho1/c
    t2 = obs2.time - rho2/c

    vel1, _ = Kepler.lambert_direct(pos1, pos2, t2 - t1, gm)
    return Kepler.Cartesian(pos1, vel1, t1, gm)
end

# assumes 0-rev lambert solution
function herget_solve_with_partials(obs1, obs2, rho1, rho2, gm, c)
    pos1 = rho1*obs1.angles + obs1.position
    pos2 = rho2*obs2.angles + obs2.position

    t1 = obs1.time - rho1/c
    t2 = obs2.time - rho2/c

    vel1, vel2 = Kepler.lambert_direct(pos1, pos2, t2 - t1, gm)

    # build the jacobian
    _, _, stm21 = Kepler.propagate_stm(pos2, vel2, t1 - t2, gm)
    dx1_dp1 = normalize(rho1*obs1.angles)
    dx2_dp2 = normalize(rho2*obs2.angles)

    dv1_dp2 = stm21.dV_dX0 * dx2_dp2

    return Kepler.Cartesian(pos1, vel1, t1, gm), dx1_dp1, dv1_dp2
end

function herget_residuals(obs, obs1, obs2, rho1, rho2, gm, c)
    orbit = herget_solve(obs1, obs2, rho1, rho2, gm, c)
    resid = compute_residuals(obs, orbit, c)
    return resid
end

function herget_residuals_with_partials(obs, obs1, obs2, rho1, rho2, gm, c; del = 1e-3)
    orbit = herget_solve(obs1, obs2, rho1, rho2, gm, c)
    resid = compute_residuals(obs, orbit, c)

    o11 = herget_solve(obs1, obs2, rho1 + del, rho2, gm, c) 
    o12 = herget_solve(obs1, obs2, rho1 - del, rho2, gm, c) 
    o21 = herget_solve(obs1, obs2, rho1, rho2 + del, gm, c)
    o22 = herget_solve(obs1, obs2, rho1, rho2 - del, gm, c)

    r11 = compute_residuals(obs, o11, c)
    r12 = compute_residuals(obs, o12, c)
    r21 = compute_residuals(obs, o21, c)
    r22 = compute_residuals(obs, o22, c)

    J1 = r11 - r12
    J2 = r21 - r22
    J  = hcat(J1, J2) ./ (2del)

    return resid, -J
end