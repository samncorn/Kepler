struct EmptySink end
Base.push!(::EmptySink, ::Any) = nothing

""" Uses Levenberg-Marquardt to solve the least squares problem

fdf
    function which returns a tuple of (function, jacobian). weighting must be done by fdf (applied to residuals and jacobian)

"""
function least_squares(fdf::F, x0;
    # mu0 = 100.0,
    # mu_a = 2.0,
    # mu_b = 4.0,
    abs_tol = 1e-8,
    rel_tol = 1e-8,
    max_iter = 100,
    # uphill_tol = 1e-3,
    sink = EmptySink(),
) where {F}
    x = copy(x0)
    # inlier_mask =

    cost, Hty, HtH = least_squares_kernel(fdf, x)
    rel_err = Inf
    abs_err = Inf
    i  = 1

    # force the initial step
    # start out not trusting gauss newton (taking a newton step)
    dx_gn = pinv(HtH)*Hty
    dx_sd = (transpose(Hty)*Hty / (transpose(Hty)*HtH*Hty)) * Hty
    mu    = norm(dx_sd + (dx_gn - dx_sd)/4)
    # mu    = norm(dx_sd)

    push!(sink, (i = i, x = x, mu = mu, cost = cost, gain = (2, 0.0)))
    while (
                 i < max_iter
        && rel_err > rel_tol
        && abs_err > abs_tol
    )
        i += 1
        # LM method
        # dx = inv(HtH + mu*diagm(diag(HtH)))*Hty
        # new_cost, new_Hty, new_HtH = least_squares_kernel(fdf, x + dx; weights = weights)

        # if new_cost > (1.0 + uphill_tol)*cost
        #     mu *= mu_a
        #     continue
        # elseif new_cost > cost
        #     mu *= mu_a
        # else
        #     mu /= mu_b
        # end

        # Powell's dog-leg method
        dx_gn = pinv(HtH)*Hty                                           # Gauss-Newton step
        dx_sd = dot(Hty, Hty) / (transpose(Hty)*HtH*Hty) * Hty          # Cauchy point
        dx_dl = dx_gn - dx_sd      
        
        # alpha = norm(dx_sd)

        dx, dL, flag = if norm(dx_gn) <= mu
            (
                dx_gn, 
                cost,
                1
            )
        elseif norm(dx_sd) < mu
            #
            a = dot(dx_dl, dx_dl)
            b = 2*dot(dx_dl, dx_sd)
            c = dot(dx_sd, dx_sd) - mu^2
            k1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
            k2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
            k  = max(k1, k2) # whichever is positive
            _dx = dx_sd + k*dx_dl
            dL_sd = 2*transpose(Hty)*dx_sd + transpose(dx_sd)*HtH*dx_sd
            _dL   = (1 - 2k)*dL_sd + k*cost
            (_dx, _dL, 2)
        else
            _dx = mu*normalize(dx_sd)
            _dL = 2*transpose(Hty)*_dx + transpose(_dx)*HtH*_dx
            (_dx, _dL, 3)
        end

        new_cost, new_Hty, new_HtH = least_squares_kernel(fdf, x + dx)
        # compare steps to the prediction, step, update trust region radius
        # should tune these values
        # 
        gain = (cost - new_cost)/abs(dL)

        if gain < 0.25 # bad approximation, shrink the trust region
            mu /= 5
            # if gain < 0.01
            #     continue # so bad we need to retry with lower mu
            # end
        elseif gain > 0.75 # good approximation, expand the trust region
            mu = max(mu, 7*norm(dx))
        end

        abs_err = abs(cost - new_cost)
        rel_err = abs_err/abs(cost)

        x += dx
        cost = new_cost
        Hty  = new_Hty
        HtH  = new_HtH
        push!(sink, (i = i, x = x, mu = mu, cost = cost, gain = (flag, gain)))
    end
    return x, cost, i
end

# weights must be handled by fdf
function least_squares_kernel(fdf::F, x) where {F}
    D = length(x)
    T = eltype(x)

    J   = zero(T)
    Hty = zeros(T, D)
    HtH = zeros(T, D, D)

    for (dyi, Hi) in fdf(x)
        if !inlier
            continue
        end
        J   += dot(dyi, dyi)
        Ht   = transpose(Hi)
        Hty += Htdyi
        HtH += HtH
    end

    return J, Hty, HtH
end

function least_squares_kernel(fdf::F, x::SVector{D, T}) where {F, D, T}
    J   = zero(T)
    Hty = @SVector zeros(T, D)
    HtH = @SMatrix zeros(T, D, D)

    for (dyi, Hi) in fdf(x)
        J   += dot(dyi, dyi)
        Ht   = transpose(Hi)
        Hty += Ht*dyi
        HtH += Ht*Hi
    end

    return J, Hty, HtH
end

""" wrapper
"""
least_squares(f::F1, df::F2, x0; kwargs...) where {F1, F2} = least_squares(x -> (f(x), df(x)), x0; kwargs...)
