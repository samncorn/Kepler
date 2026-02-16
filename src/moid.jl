""" Ellipse parameterization in the values most convenient for moid computation
"""
struct Ellipse{T, V}
    e::T # eccentricity
    p::T # semilatus rectum
    x::V # unit vector in periapse direction
    y::V # unit vector perp to x and z
end


function Ellipse(orbit::Orbit)
    hvec = cross(orbit.position, orbit.velocity)
    evec = cross(orbit.velocity, hvec)/orbit.gm - normalize(orbit.position)
    p    = dot(hvec, hvec)/orbit.gm # semilatus rectum
    e    = norm(evec)               # eccentricity
    x    = evec/e
    y    = normalize(cross(hvec, evec))
    return Ellipse(e, p, x, y)
end

function Ellipse(orbit::CometaryOrbit)
    p = orbit.q*(1.0 + orbit.e)
    # a = orbit.q/(1.0 - orbit.e)
    # p = a*(1.0 - orbit.e^2)
    x = SVector{3}(1.0, 0.0, 0.0)
    y = SVector{3}(0.0, 1.0, 0.0)
    # rotate appropriately
    R = zrot(-orbit.Om)*xrot(-orbit.i)*zrot(-orbit.w)
    return Ellipse(orbit.e, p, R*x, R*y)
end

function radius(ellipse, angle)
    return ellipse.p/(1.0 + ellipse.e*cos(angle))
end

function radial_rate(ellipse, angle)
    return ellipse.p*ellipse.e*sin(angle)/(1.0 + ellipse.e*cos(angle))^2
end

norm2(x) = dot(x, x)

""" convenience function to avoid type piracy on nothing
"""
struct EmptySink end
Base.push!(::EmptySink, _) = nothing
Base.getindex(::EmptySink, idx)  = EmptySink()
Base.lastindex(::EmptySink) = EmptySink()

""" 
"""
function moid_simplex(
    ellipse1::Kepler.Ellipse{T, V}, 
    ellipse2::Kepler.Ellipse{T, V};
    tol1::T     = 1e-4,
    tol2::T     = 1e-8, 
    tol_cond::T = 1e3,
    max_iter    = 100,
    ds::T       = 0.2,
    ) where {T, V}
    # scanning pass, get initial guesses
    _, u_min, v_min = moid_scan_meridional(ellipse1, ellipse2)

    # run levenberg-marquardt on each guess, take the lowest
    moid = Inf
    vf   = Inf
    uf   = Inf
    _f   = ((u, v),) -> sum(dists_with_partials(u, ellipse1, v, ellipse2)[1] .^ 2)
    for (v, u) in zip(v_min, u_min)
        if isinf(v) || isinf(u)
            continue
        end

        # use a simplex method to refine the initial guess
        # points = initial_simplex(SVector{2}(u, v), 0.3)
        points = MVector{3, SVector{2, T}}(
            SVector{2}(u + ds, v),
            SVector{2}(u, v + ds),
            SVector{2}(u - ds, v - ds),
        )
        vals = _f.(points)

        d0 = Inf
        x  = sum(points) ./ length(points)
        d  = _f(x)
        i  = 1
        while abs(d0 - d) > tol1 && i < max_iter
            i += 1
            simplex_step!(_f, points, vals)
            d0 = d
            x  = sum(points) ./ length(points)
            d  = _f(x)
        end
        
        # now run least squares until convergence (or worse)
        i    += 1
        u, v  = x
        dx, H = dists_with_partials(u, ellipse1, v, ellipse2)
        Ht    = transpose(H)
        d     = dot(dx, dx)
        di    = Inf
        stepx = Inf
        while stepx > tol2 && abs(d - d0) > tol2 && i < max_iter
            i += 1

            du, dv = inv(Ht*H)*Ht*dx
            dx, H  = dists_with_partials(u + du, ellipse1, v + dv, ellipse2)
            Ht     = transpose(H)
            stepx  = sqrt(du^2 + dv^2)
            di     = dot(dx, dx)

            if di > d
                simplex_step!(_f, points, vals)
                x  = sum(points) ./ length(points)
                di = _f(x)
                u, v = x

                dx, H = dists_with_partials(u, ellipse1, v, ellipse2)
                Ht    = transpose(H)
            else
                v += dv
                u += du
            end
            # if cond(Ht*H) < tol_cond
            #     # reasonably well conditioned, continue normal equations solution
            #     du, dv = inv(Ht*H)*Ht*dx
            #     dx, H = dists_with_partials(u + du, ellipse1, v + dv, ellipse2)
            #     Ht    = transpose(H)
            #     stepx = sqrt(du^2 + dv^2)
            #     di    = dot(dx, dx)

            #     if di 
                
            #     v  = mod(v + dv, 2pi)
            #     u  = mod(u + du, 2pi)
            # else
            #     # take another simplex step
            #     simplex_step!(_f, points, vals)
            #     x  = sum(points) ./ length(points)
            #     di = _f(x)

            #     u, v = x
            #     v  = mod(v, 2pi)
            #     u  = mod(u, 2pi)

            #     # still need to compute matrix for next step
            #     dx, H = dists_with_partials(u, ellipse1, v, ellipse2)
            #     Ht    = transpose(H)
            # end

            d0  = d
            d   = di
        end

        if d < moid
            moid = d
            vf   = v
            uf   = u
        end
    end

    return sqrt(moid), uf, vf
end


function simplex_step!(f, points, vals; a = 1.0, y = 2.0, r = 0.5, s = 0.5)
    j  = argmax(vals)
    xc = (sum(points) - points[j]) ./ (length(points) - 1)
    # dc = f(xc)
    x1 = xc + a*(xc - points[j])
    d1 = f(x1)

    count_better = count(val -> val < d1, vals)
    count_worse  = count(val -> val >= d1, vals)

    if count_better == 0
        # expand
        xe = xc + y*(x1 - xc)
        de = f(xe)
        if de < d1
            points[j] = xe
            vals[j]   = de
        else
            points[j] = x1
            vals[j]   = d1
        end
    elseif count_better > 0 && count_worse > 1
        # reflect
        points[j] = x1
        vals[j]   = d1
    else
        points[j] = xc + r*(x1 - xc)
        vals[j]   = f(points[j])
    end
end

""" method from Wisniowsky + Rickmann 2013. scans for near minima, then iteratively refines
"""
function moid_scan(
    ellipse1::Ellipse{T, V}, 
    ellipse2::Ellipse{T, V};
    mu::T = 10.0, 
    tol::T = 1e-8, 
    max_iter = 100, 
    mu_f::T = 5.0,
    nu::T = 1e-2,
    k::T  = 0.9,
    ) where {T, V}
    # scanning pass, get initial guesses
    dmin, u_min, v_min = moid_scan_meridional(ellipse1, ellipse2)

    # run levenberg-marquardt on each guess, take the lowest
    moid = Inf
    vf   = Inf
    uf   = Inf
    for (v, u) in zip(v_min, u_min)
        if isinf(v) || isinf(u)
            continue
        end

        dx, H = dx, H = dists_with_partials(u, ellipse1, v, ellipse2)
        d   = dot(dx, dx)
        Ht  = transpose(H)

        i = 0
        d0 = 2d
        stepx = Inf
        # k2 = k
        while abs(d - d0) > tol && stepx > tol && i < max_iter
            i += 1
            (du, dv) = pinv(Ht*H + mu*I)*Ht*dx

            dx, H = dists_with_partials(u + du, ellipse1, v + dv, ellipse2)
            stepx = abs(norm(dx))
            di  = dot(dx, dx)
            Ht  = transpose(H)

            if di >= d
                mu *= mu_f
            else
                mu /= mu_f
                nu /= mu_f
            end

            v  = mod(v + dv, 2pi)
            u  = mod(u + du, 2pi)
            d0  = d
            d   = di
        end

        if d < moid
            moid = d
            vf  = v
            uf  = u
        end
    end

    return sqrt(moid), uf, vf
end

# compute the distance between two point on orbits given a few of their paramters,
# along with the partial derivatives w.r.t true anomaly
function dists_with_partials(u, ellipse1, v, ellipse2)
    r1 = radius(ellipse1, u)
    r2 = radius(ellipse2, v)

    pos1 = r1*(cos(u)*ellipse1.x + sin(u)*ellipse1.y)
    pos2 = r2*(cos(v)*ellipse2.x + sin(v)*ellipse2.y)

    dr1 = radial_rate(ellipse1, u)
    dr2 = radial_rate(ellipse2, v)

    dx   = pos1 - pos2
    H    = hcat(
          dr1*pos1/r1 + r1*(-sin(u)*ellipse1.x + cos(u)*ellipse1.y),
        -(dr2*pos2/r2 + r2*(-sin(v)*ellipse2.x + cos(v)*ellipse2.y))
    )
    return dx, -H
end

function moid_scan_meridional(ellipse1, ellipse2)
    v_min = MVector{4}(Inf, Inf, Inf, Inf)
    u_min = MVector{4}(Inf, Inf, Inf, Inf)
    dmin  = MVector{4}(Inf, Inf, Inf, Inf)

    i     = 1
    D1    = NaN
    D2    = NaN
    v1   = NaN
    v2   = NaN
    u1   = NaN
    u2   = NaN
    for v in 0:0.12:(2pi + 0.12)
        # find the position of o2
        r2   = radius(ellipse2, v)
        pos2 = r2*(cos(v)*ellipse2.x + sin(v)*ellipse2.y)
        
        # get out of plane componet w.r.t o1
        x2  = dot(pos2, ellipse1.x)
        y2  = dot(pos2, ellipse1.y)
        z22 = r2^2 - x2^2 - y2^2

        # compute o1 radial distance for the same longitude
        rho = sqrt(x2^2 + y2^2)
        r1  = ellipse1.p/(1.0 + ellipse1.e*x2/rho)

        # compute meriodional distance
        D = sqrt(z22 + (rho - r1)^2)
        u = mod(atan(y2, x2), 2pi)

        # check for local min
        if D1 > D2 < D
            v_min[i] = v2
            u_min[i] = u2
            dmin[i]  = D2
            i += 1
        end

        D1 = D2
        D2 = D
        v1 = v2
        v2 = v
        u1 = u2
        u2 = u
    end

    if count(x -> !isinf(x), dmin) <= 1
        v_min = MVector{4}(0.0, pi/2, pi, 3pi/2)
        for (i, v) in enumerate(v_min)
            # find the position of o2
            r2   = radius(ellipse2, v)
            pos2 = r2*(cos(v)*ellipse2.x + sin(v)*ellipse2.y)
            
            # get out of plane componet w.r.t o1
            x2  = dot(pos2, ellipse1.x)
            y2  = dot(pos2, ellipse1.y)
            u_min[i] = mod(atan(y2, x2), 2pi)

            z22 = r2^2 - x2^2 - y2^2

            # compute o1 radial distance for the same longitude
            rho = sqrt(x2^2 + y2^2)
            r1  = ellipse1.p/(1.0 + ellipse1.e*x2/rho)

            # compute meriodional distance
            D = sqrt(z22 + (rho - r1)^2)
            dmin[i] = D
        end
    end

    return dmin, SVector{4}(u_min), SVector{4}(v_min)
end