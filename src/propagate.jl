function stumpff(z)
    c0, c1, c2, c3 = if z > 0
        sin2 = sin(0.5sqrt(z))
        cos2 = cos(0.5sqrt(z))
        c1 = 2sin2*cos2/sqrt(z)
        c2 = 2sin2^2/z
        c3 = (1 - c1)/z
        c0 = 1 - z*c2
        c0, c1, c2, c3
    elseif z < 0
        sinh2 = sinh(0.5sqrt(-z))
        cosh2 = cosh(0.5sqrt(-z))
        c1 = 2sinh2*cosh2/sqrt(-z)
        c2 = -2sinh2^2/z
        c3 = (1 - c1)/z
        c0 = 1 - z*c2
        c0, c1, c2, c3
    else
        c1 = 1.0
        c2 = 1/2
        c3 = 1/6
        c0 = 1.0
        c0, c1, c2, c3
    end

    return c0, c1, c2, c3
end