""" backwards compatability wrapper
"""
function solve(pos0, vel0, dt, gm; max_iter = 20)
    return Kepler.propagate(pos0, vel0, dt, gm; max_iter = max_iter)
end

""" 
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
function propagate(pos, vel, dt, gm; max_iter = 20)
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

    _, c1, c2, c3 = stumpff(a*x^2)
    y = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
    
    while sign(y) == sign(y_br)
        # take forward steps, maintaining the size of our interval
        # we could use some newton steps here, and would probably work better
        # but I'm lazy, and this is exceedingly likely to require at most 1 step
        x_br = x
        y_br = y
        x    *= 2
        _, c1, c2, c3 = stumpff(a*x^2)
        y    = x*c1 + dr0*c2*x^2 + c3*x^3 - dt0
    end

    # we now have a bracket, and can find the root.
    # chain a couple of closures together to get the right args to the right places
    # TODO: implement a bracketed root finder using higher derivatives
    bracket = (x_br, x)
    # println("bracket = $bracket")
    _f1 = _x -> (_x, stumpff(a*_x^2)...)
    _f2 = ((_x, _, _c1, _c2, _c3),) -> _x*_c1 + dr0*_c2*_x^2 + _c3*_x^3
    x = find_zero(_x -> _f2(_f1(_x)) - dt0, bracket, A42())

    _, c1, c2, c3 = stumpff(a*x^2)
    f    = 1.0 - c2*x^2
    g    = dt0 - c3*x^3
    posf = f*pos0 + g*vel0

    rf   = norm(posf)
    dg   = 1.0 - c2/rf*x^2
    df   = x*(a*c3*x^2 - 1.0)/rf
    velf = df*pos0 + dg*vel0

    # re-dimensionalize our inputs
    return posf*DU, velf*DU/TU
end