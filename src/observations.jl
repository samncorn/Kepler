struct angles_observation{T}
    time::T
    position::SVector{3, T}
    angles::SVector{3, T}
end

function apparent_position(target_traj, obs_pos, t, c; lt_tol = 1e-15, max_iter = 20)
    lt = 0.0
    dt = Inf
    i  = 0
    while dt > lt_tol && i < max_iter
        i += 1
        # println(i)
        # println(" ti = $(t - lt) ")
        lti = try
            pos = target_traj(t - lt)
            if any(isnan.(pos))
                throw((dt = t - lt))
            end
            norm(pos - obs_pos) / c
        catch _
            throw((dt = t-lt,))
        end
        dt  = abs(lt - lti)
        lt = lti
    end

    # println("evaluating light time corrected state")
    pos_app = try
        target_traj(t - lt) - obs_pos
    catch _
        throw((dt = t-lt,))
    end

    return pos_app, t - lt, i
end

function state_to_angles(x0::Kepler.Cartesian, t_obs, obs_pos, c)
    traj = t -> Kepler.propagate(x0, t)[1]
    app, _, _  = apparent_position(traj, obs_pos, t_obs, c; max_iter = 3)
    r    = norm(app)
    xhat = app/r
    return xhat
end

function state_to_angles_with_partials(x0::Kepler.Cartesian, t_obs, obs_pos, c)
    traj = t -> Kepler.propagate(x0, t).position
    app, lt, _  = apparent_position(traj, obs_pos, t_obs, c; max_iter = 3)
    r    = norm(app)
    xhat = app/r

    _, stm = Kepler.propagate_stm(x0, lt)
    dxdx = stm.dX_dX0
    dxdv = stm.dX_dV0
    
    dang_dx = (I3 - xhat*transpose(xhat))/r

    return xhat, dang_dx*dxdx, dang_dx*dxdv
end

function compute_residuals(obs, x0::Kepler.Cartesian, c)
    ang = state_to_angles(x0, obs.time, obs.position, c)
    # convert to long-lat
    z = SVector{3}(0.0, 0.0, 1.0)
    # tangent basis
    a = normalize(cross(z, obs.angles))
    d = normalize(cross(obs.angles, a))
    return -SVector{2}(dot(ang, a), dot(ang, d))
end

function compute_residuals_with_partials(obs, x0::Kepler.Cartesian, c)
    ang, dang_dx0, dang_dv0 = state_to_angles_with_partials(x0, obs.time, obs.position, c)
    # convert to long-lat + long-lat partials
    z = SVector{3}(0.0, 0.0, 1.0)
    # tangent basis
    a = normalize(cross(z, obs.angles))
    d = normalize(cross(obs.angles, a))
    # partials
    J = transpose(hcat(a, d))
    return -SVector{2}(dot(ang, a), dot(ang, d)), -J*dang_dx0, -J*dang_dv0
end