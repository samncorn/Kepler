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
function herget_kernel(rho12, observations, obs1, obs2, gm, c)
    orbit, dx1_dp1, dv1_dp2 = herget_solve(obs1, obs2, rho12..., gm, c)
    return Iterators.map(o -> herget_residuals_with_partials(o, orbit, dx1_dp1, dv1_dp2, c), observations)
end

function herget_residuals_with_partials(obs, orbit, dx1_dp1, dv1_dp2, c)
    resid, J_x, J_v = compute_residuals_with_partials(obs, orbit, c)
    return resid, J_x*dx1_dp1 + J_v*dv1_dp2
end

# assumes 0-rev lambert solution
function herget_solve_with_partials(obs1, obs2, rho1, rho2, gm, c)
    pos1 = rho1*obs1.angles + obs1.position
    pos2 = rho2*obs2.angles + obs2.position

    t1 = obs1.time - rho1/c
    t2 = obs2.time - rho2/c

    vel1, vel2 = lambert_direct(pos1, pos2, t2 - t1, gm)

    # build the jacobian
    _, _, stm21 = propagate_stm(pos2, vel2, t1 - t2, gm)
    dx1_dp1 = normalize(pos1)
    dx2_dp2 = normalize(pos2)

    dv1_dp2 = stm21.dV_dX0 * dx2_dp2

    return Kepler.Cartesian(pos1, vel1, t1, gm), dx1_dp1, dv1_dp2
end

# assumes 0-rev lambert solution
function herget_solve(obs1, obs2, rho1, rho2, gm, c)
    pos1 = rho1*obs1.angles + obs1.position
    pos2 = rho2*obs2.angles + obs2.position

    t1 = obs1.time - rho1/c
    t2 = obs2.time - rho2/c

    vel1, _ = lambert_direct(pos1, pos2, t2 - t1, gm)
    return Kepler.Cartesian(pos1, vel1, t1, gm), dx1_dp1, dv1_dp2
end