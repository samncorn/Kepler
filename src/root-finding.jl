""" f must return the appropriate amount of derivatives for next_term
"""
function root_solve(next_term, f, x0, tol, max_iter)
    i = 0
    x = x0
    err = Inf
    while err > tol && i <= max_iter
        dx = next_term(f, x)
        x -= dx
        err = abs(dx)
        i += 1
    end
    return x, i
end

function _newton2(f, x)
    y, dy = f(x)
    dx = y/dy
    return dx
end

function _newton3(f, x)
    y, dy, ddy = f(x)
    dx = y/dy
    dx2 = y/(dy - dx*ddy/2)
    return dx2
end

function _newton4(f, x)
    y, dy, ddy, dddy = f(x)
    dx = y/dy
    dx2 = y/(dy - dx*ddy/2)
    dx3 = y/(dy - dx2*ddy/2 + dx2^2 * dddy/6)
    return dx3
end

function _laguerre(f, x, n)
    y, dy, ddy = f(x)
    return n*y / (dy + sign(dy)*sqrt(abs((n-1)^2 * (dy^2 - n*(n-1)*y*ddy))))
end

newton2(f, x0, tol, max_iter) = root_solve(_newton2, f, x0, tol, max_iter)
newton3(f, x0, tol, max_iter) = root_solve(_newton3, f, x0, tol, max_iter)
newton4(f, x0, tol, max_iter) = root_solve(_newton4, f, x0, tol, max_iter)
# laguerre(f, x0, tol, max_iter; n = 5) = root_solve(x -> _laguerre(x[1], x[2], 5), f, x0, tol, max_iter)
function laguerre(f, x0, tol, max_iter; n = 5)
    i = 0
    x = x0
    err = Inf
    while err > tol && i <= max_iter
        dx = _laguerre(f, x, n)
        x -= dx
        err = abs(dx)
        i += 1
    end
    return x, i
end