# function regula_falsi(f, x0, tol, max_iter)
#     i = 0
#     x = x0
#     err = Inf
#     while err > tol && i <= max_iter
#         y = f(x)
#         err = abs(y - x)
#         x = y
#         i += 1
#     end
#     return x, i
# end

function newton2(f, x0, tol, max_iter)
    i = 0
    x = x0
    err = Inf
    while err > tol && i <= max_iter
        y, dy = f(x)
        dx = y/dy
        err = abs(dx)
        x -= dx
        i += 1
    end
    return x, i
end

function newton3(f, x0, tol, max_iter)
    i = 0 
    x = x0
    err = Inf
    while err > tol && i <= max_iter
        y, dy, ddy = f(x)
        dx = y/dy
        dx2 = y/(dy - dx*ddy/2)
        x -= dx2
        err = abs(dx2)
        i += 1
    end
    return x, i
end

function newton4(f, x0, tol, max_iter)
    i = 0
    x = x0
    err = Inf
    while err > tol && i <= max_iter
        y, dy, ddy, dddy = f(x)
        dx = y/dy
        dx2 = y/(dy - dx*ddy/2)
        dx3 = y/(dy - dx2*ddy/2 + dx2^2 * dddy/6)
        err = abs(dx3)
        x -= dx3
        i += 1
    end
    return x, i
end