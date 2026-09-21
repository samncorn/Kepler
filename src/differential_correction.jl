""" Fits a keplerian orbit to a set of optical measurements.

Arguments
=========
observations
    iterator of measurements. Each element must have the fields
        'angles'    ... cartesian unit vector of ra, dec
        'time'      ...
        'position'  ... observer position
init_state
    Cartesian state describing the intitial guess of the orbit. Must have the fields
        'position'
        'velocity'
        'epoch'
        'gm'

kwargs
======
c (defualt = Inf)
    speed of light, for light time correction.
kwargs
    least sqaures key word arguments (see Kepler.least_squares)

"""
function fit_kepler(
    observations,
    init_state;
    c = Inf,
    kwargs...
)
    # TODO
    # make generic over dynamical model
        # the initial state needs a method that propagates to meas times (w/ partial derivatives)
        #
    # make generic over measurement type (optical, radar, etc)
        # Don't actually need to implement anything other than optical for now, just iron out the interface

    x0 = reduce(vcat, (init_state.position, init_state.velocity))
    xf, chi2, i = Kepler.least_squares(x -> Iterators.map(
            o -> compute_weighted_residuals_with_partials(
                o,
                Kepler.Cartesian(
                    SVector{3}(x[1], x[2], x[3]),
                    SVector{3}(x[4], x[5], x[6]),
                    init_state.epoch,
                    init_state.gm
                ),
                c
            ),
            observations
        ),
        x0;
        kwargs...
    )

    cart = Kepler.Cartesian(
        SVector{3}(xf[1], xf[2], xf[3]), 
        SVector{3}(xf[4], xf[5], xf[6]), 
        init_state.epoch, 
        init_state.gm
    )

    # TODO: return post-fit covariance matrix?
    return cart, chi2, i
end

function fit_kepler_progressive_rejection(observations, initial_orbit; 
    c = Inf, 
    min_obs = 6, 
    max_rej = 1, 
    X2_rej = 8.0, 
    X2_rec = 7.0, 
    alpha = 0.25,
    # beta = 0.1,
    LS_kwargs...
)
    inliers         = trues(length(observations))
    n_rem           = length(observations)
    n_rej           = 0
    X2_max          = 0.0
    X2_max_running  = Inf
    x0 = reduce(vcat, (initial_orbit.position, initial_orbit.velocity))

    j = 0

    function get_resid(o, x)
        return Kepler.compute_weighted_residuals_with_partials(
            o,
            Kepler.Cartesian(
                SVector{3}(x[1], x[2], x[3]),
                SVector{3}(x[4], x[5], x[6]),
                initial_orbit.epoch,
                initial_orbit.gm
            ),
            c
        )
    end

    function get_resid_masking(o, x, i)
        dy, H = Kepler.compute_weighted_residuals_with_partials(
            o,
            Kepler.Cartesian(
                SVector{3}(x[1], x[2], x[3]),
                SVector{3}(x[4], x[5], x[6]),
                initial_orbit.epoch,
                initial_orbit.gm
            ),
            c
        )

        X2 = dot(dy, dy)

        # X2_max_running = max(X2, X2_max_running)
        
        if inliers[i]
            # check if we should exclude the obs
            if X2 > max(X2_rej, alpha*X2_max) && n_rej < max_rej && n_rem > min_obs
                n_rej += 1
                n_rem -= 1
                inliers[i] = false
                # println("rejecting observation $i with X2 = $X2 > $(max(X2_rej, alpha*X2_max)) (n_rej = $n_rej, n_rem = $n_rem)")
            end
        else
            # check whether we want to recover
            if X2 < X2_rec
                inliers[i] = true
                n_rem += 1
                # println("recovering observation $i (n_rem = $n_rem)")
            end
        end

        if inliers[i]
            return (dy, H)
        else
            # if at the end of the checks we are excluding this observation,
            # zero it out
            return(dy - dy, H - H)
        end
    end

    function get_resids(x)
        # reset the trackers 
        n_rej -= n_rej
        # X2_max          = X2_max_running
        # X2_max_running  = 0.0
        j += 1
        X2_max = 0.0
        
        for (i, (dy, _)) in enumerate(Iterators.map(o -> get_resid(o, x), observations))
            if inliers[i]
                X2 = dot(dy, dy)
                X2_max = max(X2_max, X2)
            end
        end
        # println("iteration $j, max X2 = $X2_max")
        
        return Iterators.map(((i, o),) -> get_resid_masking(o, x, i), enumerate(observations))
    end


    xf, X2, i = Kepler.least_squares(get_resids, x0; LS_kwargs...)
    cart = Kepler.Cartesian(SVector{3}(xf[1], xf[2], xf[3]), SVector{3}(xf[4], xf[5], xf[6]), initial_orbit.epoch, initial_orbit.gm)
    return cart, X2, i, inliers
end
