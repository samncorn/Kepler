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
    # fdf = x -> Iterators.map(((t, o, ang),) -> kepler_resids_with_partials_analytic(x, t-t0, o, ang, gm, c), zip(times, obs, angles))
    x0 = reduce(vcat, (init_state.position, init_state.velocity))
    xf, _, _ = Kepler.least_squares(x -> Iterators.map(
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

    # TODO: return post-fit covariance matrix?
    return Kepler.Cartesian(
        SVector{3}(xf[1], xf[2], xf[3]),
        SVector{3}(xf[4], xf[5], xf[6]),
        init_state.epoch,
        init_state.gm
    )
    # levenberg_marquardt(fdf, x0, weights; mu = mu, tol_cost_rel = tol, tol_cost_abs = tol, tol_param_rel = tol, tol_param_abs = tol, max_iter = max_iter, sink = sink)
end

mutable struct inlier_cache{T}
    residuals::Vector{T}
    inlier_mask::BitArray
    iteration_rejected::Int
end

# inlier_cache(observations) = inlier_cache(zeros(), trues(length(observations)), 0, length(observations))

"""
"""
# TODO: Hook in callbacks for setting rejection flags and checking mutually exclusive observations
#       Residual comparisons require sparse pairwise combinations to be checked. So it is probably
#       better to save residuals to a preallocated array on each iteration. This should have negligible 
#       performance impact, since the write time should be negligible compared to the proapgation time,
#       even for 2-body solves.
#       Also, callbacks to terminate on 
function fit_analytic_reject_outliers(
    observations,
    init_state;
    c = Inf,
    # inliers = inlier_cache(observations),
    inlier_mask = trues(length(observations)),
    resid_cut = 5.0, iteration_rejection_max = 10, min_inliers = 6,
    kwargs...
)
    # # initialize the inlier set by computing the initial residuals
    # x0 = reduce(vcat, (init_state.position, init_state.velocity))
    # xf, _, _ = Kepler.least_squares(o -> Iterators.map(
    #             ((i, x),) -> set_inlier!(inlier_mask, i, compute_weighted_residuals_with_partials(
    #                     o,
    #                     Kepler.Cartesian(
    #                         SVector{3}[x[1], x[2], x[3]],
    #                         SVector{3}[x[4], x[5], x[6]],
    #                         init_state.epoch,
    #                         init_state.gm
    #                     ),
    #                     c
    #                 );
    #                 resid_cut = resid_cut,
    #                 iteration_rejection_max = iteration_rejection_max,
    #                 min_inliers = min_inliers
    #             ),
    #             enumerate(observations)
    #         ),
    #     ),
    #     x0;
    #     kwargs...
    # )

    # # TODO: return post-fit covariance matrix?
    # return Kepler.Cartesian(
    #     SVector{3}[xf[1], xf[2], xf[3]],
    #     SVector{3}[xf[4], xf[5], xf[6]],
    #     init_state.epoch,
    #     init_state.gm
    # )
end

function set_inlier!(inliers, i, resid_with_partial; resid_cut = 5.0, iteration_rejection_max = 10, min_inliers = 6)
    if i == 1
        inliers.iteration_rejected == 0
    end

    if i == length(inliers.inlier_mask)
        # run mutual exclusivity callback

    end

    resid = resid_with_partial[1]
    if dot(resid, resid) > resid_cut && sum(inliers.inlier_mask) > min_inliers
        inliers.inlier_mask[i] = false
        inliers.iteration_rejected += 1
    end
    return resid_with_partial
end
