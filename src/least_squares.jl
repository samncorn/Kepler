struct EmptySink end
Base.push!(::EmptySink, ::Any) = nothing

""" Uses Levenberg-Marquardt to solve the least squares problem

fdf is a function which returns a tuple of (function, jacobian)
"""
function least_squares(fdf, x0; 
    weights = Iterators.cycle((I,)), 
    mu0 = 100.0, 
    mu_a = 2.0, 
    mu_b = 4.0, 
    abs_tol = 1e-8, 
    rel_tol = 1e-8, 
    max_iter = 100, 
    uphill_tol = 1e-3,
    sink = EmptySink(),
)
    x = copy(x0)

    cost, Hty, HtH = least_squares_kernel(fdf, x; weights = weights)
    rel_err = Inf
    abs_err = Inf
    i  = 0
    mu = mu0

    push!(sink, (i = i, x = x, mu = mu, cost = cost))
    while (
                 i < max_iter
        && rel_err > rel_tol
        && abs_err > abs_tol
    )
        i += 1
        dx = inv(HtH + mu*diagm(diag(HtH)))*Hty
        
        new_cost, new_Hty, new_HtH = least_squares_kernel(fdf, x + dx; weights = weights)

        if new_cost > (1.0 + uphill_tol)*cost
            mu *= mu_a
            continue
        elseif new_cost > cost
            mu *= mu_a
        else
            mu /= mu_b
        end

        abs_err = abs(cost - new_cost)
        rel_err = abs_err/abs(cost)

        x += dx
        cost = new_cost
        Hty  = new_Hty
        HtH  = new_HtH
        push!(sink, (i = i, x = x, mu = mu, cost = cost))
    end
    return x, cost, i
end

function least_squares_kernel(fdf, x; weights = Iterators.cycle((I,)))
    D = length(x)
    T = eltype(x)

    J   = zero(T)
    Hty = zeros(T, D)
    HtH = zeros(T, D, D)

    for ((dyi, Hi), wi) in zip(fdf(x), weights)
        J   += dot(dyi, wi*dyi)
        Ht   = transpose(Hi)
        Hty += Ht*wi*dyi
        HtH += Ht*wi*H
    end

    return J, Hty, HtH
end

function least_squares_kernel(fdf, x::SVector{D, T}; weights = Iterators.cycle((I,))) where {D, T}
    J   = zero(T)
    Hty = @SVector zeros(T, D)
    HtH = @SMatrix zeros(T, D, D)

    for ((dyi, Hi), wi) in zip(fdf(x), weights)
        J   += dot(dyi, wi*dyi)
        Ht   = transpose(Hi)
        Hty += Ht*wi*dyi
        HtH += Ht*wi*Hi
    end

    return J, Hty, HtH
end

""" wrappers
"""
least_squares(f, df, x0; kwargs...) = least_squares(x -> (f(x), df(x)), x0; kwargs...)