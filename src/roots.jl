function chandrupatla(f::F, bracket) where {F}
    t = 0.5
    x2, x3 = bracket

    y2 = f(x2)
    y3 = f(x3)

    if !(sign(y1[1]*y2[1]) <= 0)
        throw("Not a bracket")
    end

    x1 = x2 + t*(x3 - x2)
    i = 0
    dx1 = abs(x3 - x2)
    dx2 = 2dx1
    dx3 = 2dx2
    while x1 != x2 && x1 != x3
        if i == 100
            throw("max iters")
        end
        i += 1

        y1 = f(x1)

        if sign(y1) == sign(y2)
            (x2, x3, y2, y3) = (x3, x2, y3, y2)
        end

        phi = (y1 - y2)/(y3 - y2)
        xi  = (x1 - x2)/(x3 - x2)

        if 4dx1 < dx3 && phi^2 < xi < 1 - (1 - phi)^2
            a = (x3 - x1)/(x2 - x1)
            t = (y1/(y1 - y2))*(y3/(y3 - y2)) - a*(y1/(y3 - y1))*(y2/(y2 - y3))
            # t *= 0.95
        else 
            t = 0.5
        end
        
        # if abs(y1) < 1e-13
        #     # safeguard against too good of a guess
        #     # t = min(0.5, 4*t)
        #     t = 0.25
        # end
        
        (x3, y3, dx3) = (x2, y2, dx2)
        (x2, y2, dx2) = (x1, y1, dx1)
        x1  = x2 + t*(x3 - x2)
        dx1 = abs(t*(x3 - x2))

    end
    return x1
end

function newton_bracket(f, bracket)
    x1, x2 = bracket
    y1 = f(x1)
    y2 = f(x2)
    
    if y2[1] > y1[1]
        (x2, x1) = (x1, x2)
        (y2, y1) = (y1, y2)
    end

    x   = (x2*y2[2] - y2[1])/y2[2]
    dx1 = abs(x - x2)
    dx2 = abs(x2 - x1)
    dx3 = 2dx2
    i = 0
    while x != x1 && x != x2
        if i == 400
            throw("max iters")
        end
        i += 1
        y = f(x)
        if sign(y[1]) == sign(y2[1])
            (x1, x2, y1, y2) = (x2, x1, y2, y1)
        end

        _x = if 4dx1 < dx3
            (x*y[2] - y[1])/y[2]
        else
            x + (x2 - x)/2 
        end

        if min(x, x2) > _x || max(x, x2) < _x
            _x =  x + (x2 - x)/2
        end

        (x1, y1) = (x2, y1)
        (x2, y2) = (x, y)
        x = _x
        dx3, dx2, dx1 = (dx2, dx1, abs(x2 - x))
    end
    return x
end

# function laguerre_bracket(f, bracket)

# end


# "experimental"
function inverse_cubic_2point(f, bracket)
    max_iter = 100 # hard code this
    x1, x2 = bracket

    p1 = f(x1)
    p2 = f(x2)
    
    if !(sign(p1[1]*p2[1]) <= 0)
        throw("Not a bracket")
    end

    # println()
    # println("START ")
    # println()
    # println("x = [$x1, $x2]")
    # println("f = [$(p1[1]), $(p2[1])]")
    
    # force a bisection to catch too good of a guess
    x = x1 + (x2 - x1)/2
    # for i in 1:max_iter

    dx1 = abs(x - x1)
    dx2 = 2dx1
    dx3 = 2dx2
    i = 0
    while x != x1 && x != x2
        if i == max_iter
            throw("failed to converge in less than $(max_iter) iterations") 
        end
        i += 1
        
        p = f(x)

        if p[1] == 0
            return x
        end
        
        if !(sign(p[1]) == sign(p2[1]))
            # swap bracket ordering so that [x, x1] is a bracket
            (x1, x2, p1, p2) = (x2, x1, p2, p1)
        end

        phi = (p[1] - p1[1])/(p2[1] - p1[1])
        xi  = (x - x1)/(x2 - x1)

        # println()
        # println("x = [$x, $x1]")
        # println("f = [$(p[1]), $(p1[1])]")

        _x = if 4dx1 > dx3 !(phi^2 < xi && xi < (1 - (1-phi)^2))
            _x = x + (x1 - x)/2
            # println("bisect x")
            _x
        else
            _x = inverse_cubic_2point_step((x, p...), (x1, p1...))
            if _x < min(x, x1) || _x > max(x, x1)
                _x = x + (x1 - x)/2
            end
            # println("interpolate")
            _x
        end

        dx1, dx2, dx3 = (abs(_x - x), dx1, dx2)
        (x2, p2) = (x1, p1)
        (x1, p1) = (x, p)
        x = _x
    end
    return x
end

function inverse_cubic_2point_step(p1, p2)
    q = p1[2]/p2[2]

    a1 = (1 - 3q)/(q - 1)^3
    a2 = 1 + a1
    b1 = q/(q - 1)^2
    b2 = q*b1

    return a2*p2[1] - a1*p1[1] - p2[2]*(b2/p2[3] + b1/p1[3])
end
