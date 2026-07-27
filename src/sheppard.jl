# TODO:
# perf analysis of
# initial guess
# argument folding vs trig expressions
# continued fractions?

# uses the Sheppard 1984 method of continued fractions for shorter timespans
# longer timespans, convert to the trig expression for the sheppard independant variable
# Note that for evaluating universal expressions, b for each orbit is constant so the branch predictor (or even inlining) should handle it well
function propagate_sheppard(pos, vel, dt, gm)
    # non-dimensional_parameters
    # there are many ways to choose length and time scales
    # Fukushima's is undesirable for near collision orbits, which are interesting for flybys
    # energy or momentum based scales can lead to division by zero or parabolic or collisoin orbits
    # we first choose our scales to be canonical, so that gm => 1, which gives us TU = sqrt(DU^3 / gm)
    # the length scale is then arbitrary, with the restriction that it be well behaved in all cases.
    # the most obvious choice is then DU = r0, though this means that the scale is not constant across the orbit
    # but I suspect this to not be an issue
    r0 = norm(pos)
    ts = sqrt(r0^3 / abs(gm))
    T = dt / ts # non-dim time of flight
    b = 2.0 - dot(vel, vel)*r0/gm   # non-dim negative energy
    k = dot(vel, pos) / sqrt(r0*gm) # non-dim r . v

    # non-dim kepler is then T = U1(s, g) + k*U2(s, g) + U3(s, g)
    # w/ dT/ds = r/r0
    # we need to solve for s
    # instead of directly iterating on s, we can do another change of variables to u = U1(s/4, g)/U0(s/4, g) = sqrt(g)tan(4(sqrt(g)*s))
    # Sheppard gave expressions for U0(s, g), U1(s, g), U2(s, g) from u, so that only U3 needs to be found
    # for near parabolic and short times, the continued fraction evaluation of the hypergeometric function is most efficient
    # but for large times the recursion relation is more efficient and well defined, requiring a single atan call per iteration.

    # for periodic orbits, we are restricted to a single revolution
    P  = 2pi/sqrt(abs(b)^3)
    dT = T + P/2 - 2k/b
    n  = 0
    if b > 0
        n = floor(Int, dT/P)
        T -= n*P
    end

    u = solve_kepler_sheppard(T, b, k)

    # we can now compute the final state as a linear combination of the initial states
    _, U1, U2, _ = universal_sheppard(u, b)

    f = -U2
    g = (U1 + k*U2)*ts
    posf = f*pos + g*vel + pos

    rf = norm(posf)/r0
    df = -U1/rf/ts
    dg = -U2/rf
    velf = df*pos + dg*vel + vel

    return posf, velf
end

function propagate_sheppard_with_partials(pos, vel, dt, gm)
    throw("doesn't work yet")
    # non-dimensional_parameters
    # there are many ways to choose length and time scales
    # Fukushima's is undesirable for near collision orbits, which are interesting for flybys
    # energy or momentum based scales can lead to division by zero or parabolic or collisoin orbits
    # we first choose our scales to be canonical, so that gm => 1, which gives us TU = sqrt(DU^3 / gm)
    # the length scale is then arbitrary, with the restriction that it be well behaved in all cases.
    # the most obvious choice is then DU = r0, though this means that the scale is not constant across the orbit
    # but I suspect this to not be an issue
    r0 = norm(pos)
    ts = sqrt(r0^3 / abs(gm))
    T = dt / ts # non-dim time of flight
    b = 2.0 - dot(vel, vel)*r0/gm   # non-dim negative energy
    k = dot(vel, pos) / sqrt(r0*gm) # non-dim r . v

    # non-dim kepler is then T = U1(s, g) + k*U2(s, g) + U3(s, g)
    # w/ dT/ds = r/r0
    # we need to solve for s
    # instead of directly iterating on s, we can do another change of variables to u = U1(s/4, g)/U0(s/4, g) = sqrt(g)tan(4(sqrt(g)*s))
    # Sheppard gave expressions for U0(s, g), U1(s, g), U2(s, g) from u, so that only U3 needs to be found
    # for near parabolic and short times, the continued fraction evaluation of the hypergeometric function is most efficient
    # but for large times the recursion relation is more efficient and well defined, requiring a single atan call per iteration.

    # for periodic orbits, we are restricted to a single revolution
    P  = 2pi/sqrt(abs(b)^3)
    dT = T + P/2 - 2k/b
    n  = 0
    if b > 0
        n = floor(Int, dT/P)
        T -= n*P
    end
    dU  = b > 0 ? n*P / b : zero(T)

    u = solve_kepler_sheppard(T, b, k)

    # we can now compute the final state as a linear combination of the initial states
    _, U1, U2, U3, U = universal_sheppard_stm(u, b)
    U += dU

    f = -U2
    g = (U1 + k*U2)*ts
    posf = f*pos + g*vel + pos

    rf = norm(posf)/r0
    df = -U1/rf/ts
    dg = -U2/rf
    velf = df*pos + dg*vel + vel

    # C    = 3U5 - x*U4 - dt*U2
    C = U2*U3 - 3*U - (T + n*P)*U2
    dx_dx0 = stm_pos_pos0_normalized(pos/r0, posf/r0, vel/r0*ts, velf/r0*ts, rf/r0, f, C)        # units of DU/DU = 1
    dx_dv0 = stm_pos_vel0_normalized(pos/r0, posf/r0, vel/r0*ts, velf/r0*ts, f, g/ts, C)         # Units of TU
    dv_dx0 = stm_vel_pos0_normalized(pos/r0, posf/r0, vel/r0*ts, velf/r0*ts, rf/r0, df*ts/r0, C) # Units of 1/TU
    dv_dv0 = stm_vel_vel0_normalized(pos/r0, posf/r0, vel/r0*ts, velf/r0*ts, rf/r0, f, dg*r0, C) # Units of DU/TU*TU/DU = 1

    # return posf*DU, velf*DU/TU, dxdx, dxdv*TU, dvdx/TU, dvdv
    return posf, velf, KeplerSTM(dx_dx0, dx_dv0*ts, dv_dx0/ts, dv_dv0)
end

function solve_kepler_sheppard(T, b, k)
    # benchmarks indicate this method to be only marginally faster
    # thought this is on an apple silicon m2, so it may be down to hardware trig functions
    
    # initial guess
    # probably have a few options
    # if our timestep is appropriate, there is one single solution and T(u) is monotonically increasing
    # so even the most basic of methods will eventually converge as long as we're careful to actually establish a bracket
    # note that dt/du = (4*(1 - q)/r) and q = b*u^2 / (1 + b*u^2) < 1/2 (trust me)
    # a better guess is provided by Battin, but is still based on an expansion at u = 0, which is undesirable for parabolic and hyperbolic orbits with long timesteps
    # but long time step solutions will be found eventually anyways, since we can bracket the root, so it shouldn't be an issue
    # is there a way to better estimate the long term behavior? (check the vallado guesses) 
    # 
    u0 = 0.0
    t0, df0 = tof_nondim(u0, b, k)
    f0 = t0 - T

    # this one is effectively a single newton step
    s1 = T

    # from battin, a 4th order taylor expansion
    # unfortunately, the non polynomial behavior of the time-of-flight means this isnt great
    # s1 = sign(T)*abs(T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4))

    # attempts to use vallados guesses. The parabolic one seems to work fine, but the hyperbolic one does not
    # s1 = if b > 1e-6
    #     # sign(T)*abs(T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4))
    #     T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4)
    # else
    #     if b < 1e-6
    #         sign(T)*log((-2*b*T)/(k + sign(T)*(1 - b)/sqrt(-b)))/sqrt(-b)
    #         # sign(T)*log(abs(2b*T/(k - (1 - b)/sqrt(-b))))/sqrt(-b)
    #     else
    #         hvec = cross(pos, vel)
    #         p = dot(hvec, hvec)
    #         z = acot(3*T*sqrt(1/p^3))/2
    #         w = atan(cbrt(tan(z)))
    #         2*sqrt(p)*cot(2w)
    #     end
    # end

    u1 = if b > 0
        tan(sqrt(b)*abs(s1)/4)/sqrt(b)
    elseif b < 0
        tanh(sqrt(-b)*abs(s1)/4)/sqrt(-b)
    elseif b == 0
        abs(s1)/4
    else
        throw("invalid b")
    end
    u1 = sign(T)*min(abs(u1), 1/sqrt(abs(b))) # we can bound u just in case the guess overshoots wildly

    if isnan(u1)
        throw("NAN")
    end

    t1, df1 = tof_nondim(u1, b, k)
    f1 = t1 - T

    # we know s = 0 => u = 0, so f(u = 0) = 0
    # then we get a second value form battin's series reversion
    # if that's not a bracket, the best we can do is shift around until we have one
    # or we exhast the domain, in which case the input was invalid
    i = 0
    while sign(f0)*sign(f1) > 0 
        i += 1
        if i == 1000
            throw("max iterations")
        end
        if isfinite(f1)
            # initial guess fell short, we can shift
            (t0, f0, df0) = (t1, f1, df1)
            u1 = (u1 - u0) + u1
        elseif isfinite(f0) 
            # overshot so far that the guess overflows, have to scale back
            u1 = (u0 + u1)/2
        else
            # both are infinite, we cannot find the bracket
            throw(DomainError(dt, "The time span exceeds the provided precision for this orbit. Consider normalizing to a different scale or taking intermediate steps."))
        end
        t1, df1 = Kepler.tof_nondim(u1, b, k)
        f1 = t1 - T
    end

    # We can now use a host of root finding methods, and as long as we always shrink the bracket we're guarunteed to converge
    # but since we have derivative information and a bracket, we can do better
    # by using cubic hermite interpolations of the inverse, we have a convergence rate of ~2.73
    # by using a third point (quintic interpolation) we could get 2.91, but the evaluation is increasing complex for marginal gains
    # 2 points with 2nd derivatives may be worth it though (3.79)

    # this formulation for the step is (hopefully) more numerically stable
    # ui = u0 == u1 ? u0 : flmsm1_step((x = u0, y = f0, dy = df0), (x = u1, y = f1, dy = df1))
    ui = if u0 == u1
        u0
    elseif f1 == 0
        u1
    else
        # flmsm1_step((x = u0, y = f0, dy = df0), (x = u1, y = f1, dy = df1))
        dfi = (f1 - f0)/(u1 - u0)
        ui  = u1 - f1/dfi
    end

    i = 0
    du2 = Inf
    du1 = abs(u0 - u1)
    while ui != u0 && ui != u1
        i += 1
        if i == 1000
            throw("max iterations")
        end
        ti, dfi = tof_nondim(ui, b, k)
        fi = ti - T
        # shift appropriate end point
        if fi == 0
            break
        elseif sign(fi) == sign(f0)
            u0, f0, df0 = (ui, fi, dfi)
        elseif sign(fi) == sign(f1)
            u1, f1, df1 = (ui, fi, dfi)
        else
            throw("invalid sign detected $(fi), iteration $i, b = $b, T/P = $(T/P)")
        end
        
        # how much is the interval shrinking
        du0 = abs(u0 - u1)
        
        # compute next point
        
        # inverse cubic hermite interpolation (seems to be quite slow?)
        ui = Kepler.flmsm1_step((x = u0, y = f0, dy = df0), (x = u1, y = f1, dy = df1))

        # DISREGARD, the derivative calculationw as wrong, leading to bad performance by the interpolation
        # secant step, appears to be faster (results in less calls to bisection)
        # the step selector may be poorly conditioned
        # dfi = (f1 - f0)/(u1 - u0)
        # ui -= fi/dfi
        
        if !(min(u0, u1) <= ui <= max(u0, u1)) || du0 > du2/4 # to rescue us from arbitrailly long convergence we fall back to a bisection
            ui = (u0 + u1)/2
        end

        du2 = du1
        du1 = du0
    end

    return ui
end


""" non-dimensionalized (r0 = 1, gm = 1) Kepler's problem in the form given by Sheppard 1984. A variation on the universal form, but with a change of independant variable.
"""
function tof_nondim(u, b, k)
    U0, U1, U2, U3 = universal_sheppard(u, b)
    q  = b*u^2 / (1 + b*u^2)
    f  = U1 + k*U2 + U3
    r  = U0 + k*U1 + U2
    df = 4*(1 - q)*r
    return f, df
end

# by using the sheppard u independant variable, we need only one (but still 1) trig or inverse trig evaluation per iteration
# 
function universal_sheppard(u, b)
    q  = b*u^2 / (1 + b*u^2)
    V0 = (1 - 2q)
    V1 = 2*(1 - q)*u
    U0 = 2*V0^2 - 1
    U1 = 2*V0*V1
    U2 = 2*V1^2

    U3 = if abs(b) > 1e-4 # just need to avoid 1/g for small g. 
        s = if b > 0
            4*atan(u*sqrt(b))/sqrt(b)
        else
            4*atanh(u*sqrt(-b))/sqrt(-b)
        end
        (s - U1)/b
    else # should converge reasonable fast for small g
        # U = 16/15*V1^5*sheppard_continued_fraction(q)
        4/3*V1^3*sheppard_continued_fraction(3.0, 0.0, 1.5, q)
    end

    return U0, U1, U2, U3
end

function universal_sheppard_stm(u, b)
    q  = b*u^2 / (1 + b*u^2)
    V0 = (1 - 2q)
    V1 = 2*(1 - q)*u
    U0 = 2*V0^2 - 1
    U1 = 2*V0*V1
    U2 = 2*V1^2

    U3, U = if abs(b) > 1e-4 # just need to avoid 1/g for small g. 
        s = if b > 0
            4*atan(u*sqrt(b))/sqrt(b)
        else
            4*atanh(u*sqrt(-b))/sqrt(-b)
        end
        _U3 = (s     - U1)/b
        _U  = (_U3 - U1*U2/3)/b
        (
            _U3,
            _U
        )
    else # should converge reasonable fast for small g
        _U  = 16/15*V1^5*sheppard_continued_fraction(5.0, 0.0, 3.5, q)
        _U3 = b*U + U1*U2/3
    end

    return U0, U1, U2, U3, U
end

# function universal_sheppard_stm(u, b)
    # q = b*u^2 / (1 + b*u^2)

    # V0 = (1 - 2q)
    # V1 = 2*(1 - q)*u

    # U0 = 2*V0^2 - 1
    # U1 = 2*V0*V1
    # U2 = 1*V1^2

    # U3 = if abs(q) > 0.01 
    #     s = if b > 0 # very annoying
    #         atan(u/sqrt(b))/sqrt(b)
    #     else
    #         atanh(u/sqrt(-b))/sqrt(-b)
    #     end
    #     (1.0 - s*U1)/b
    # else # time elsewhere, about as fast as an atan call
    #     4/3*V1^2*gauss_continued_fraction(3, 0, 3/2, q)
    # end

    # return U0, U1, U2, U3
# end

function sheppard_continued_fraction(a::T, b::T, c::T, x::T) where {T} 
    k::T = 1 - 2*(a - b)
    l::T = 2*(c - 1)
    d::T = 4*c*(c - 1)
    n::T = 4*b*(c - a)
    A::T = 1.0
    B::T = 1.0
    G::T = 1.0
    G0 = NaN
    # i = 0
    while G != G0
    # for i in 1:10 # sufficient for |x| <= 1e-2 (possibly more)
        # i += 1
        G0 = G
        k *= -1
        l += 2
        d += 4*l
        n += (1 + k)*l
        A = d/(d - n*A*x)
        B *= A - 1
        G += B
    end
    return G
    # return i
end