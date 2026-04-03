function stumpff(z::T) where {T}
    tol = 1.0
    if z > tol
        sin2 = sin(sqrt(z)/2)
        cos2 = cos(sqrt(z)/2)

        c1 = 2sin2*cos2/sqrt(z)
        c2 = 2sin2^2/z
        c0 = 1 - z*c2
        c3 = (1 - c1)/z

        return c0, c1, c2, c3
    elseif z < -tol
        sin2 = sinh(sqrt(-z)/2)
        cos2 = cosh(sqrt(-z)/2)

        c1 = 2sin2*cos2/sqrt(-z)
        c2 = -2sin2^2/z
        c0 = 1 - z*c2
        c3 = (1 - c1)/z

        return c0, c1, c2, c3
    else
        c3 = stumpff_series(3, z)
        c2 = stumpff_series(2, z)
        c1 = 1 - z*c3
        c0 = 1 - z*c2
        return c0, c1, c2, c3
    end
end

function stumpff5(z)
    tol = 1.0
    if z > tol
        c0 = cos(sqrt(z))
        c1 = sin(sqrt(z))/sqrt(z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        c4 = (1/2 - c2)/z
        c5 = (1/6 - c3)/z
        return c0, c1, c2, c3, c4, c5
    elseif z < -tol
        c0 = cosh(sqrt(-z))
        c1 = sinh(sqrt(-z))/sqrt(-z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        c4 = (1/2 - c2)/z
        c5 = (1/6 - c3)/z
        return c0, c1, c2, c3, c4, c5
    else
        c5 = stumpff_series(5, z)
        c4 = stumpff_series(4, z)
        c3 = 1/6 - z*c5
        c2 = 1/2 - z*c4
        c1 = 1 - z*c3
        c0 = 1 - z*c2
        
        return c0, c1, c2, c3, c4, c5
    end
end

function stumpff_series(i, z::T) where {T}
    c = T(1/factorial(i))
    p = c
    d = 2c
    j = 1
    while c != d
        d = c
        p *= -z/((i + 2j)*(i + 2j - 1))
        c += p
        j += 1
    end
    return c
end

function universal03(b, s)
    if b > 0
        z  = sqrt(b)*s
        s2 = sin(z/2)
        c2 = cos(z/2)

        U1 = 2s2*c2/sqrt(b)
        U2 = 2s2*s2/b
        U0 = 1.0 - b*U2
        U3 = (s - U1)/b

        return U0, U1, U2, U3
    elseif b < 0
        z  = sqrt(-b)*s
        s2 = sinh(z/2)
        c2 = cosh(z/2)

        U1 = 2s2*c2/sqrt(-b)
        U2 = -2s2*s2/b
        U0 = 1.0 - b*U2
        U3 = (s - U1)/b

        return U0, U1, U2, U3
    else
        return 1.0, s, (s^2)/2, (s^3)/6
    end
end

function universal05(b, s)
    if b > 0
        z  = sqrt(b)*s
        s2 = sin(z/2)
        c2 = cos(z/2)

        U1 = 2s2*c2/sqrt(b)
        U2 = 2s2*s2/b
        U0 = 1.0 - b*U2
        U3 = (s - U1)/b
        U4 = (s^2/2.0 - U2)/b
        U5 = (s^3/6.0 - U3)/b

        return U0, U1, U2, U3, U4, U5
    elseif b < 0
        z  = sqrt(-b)*s
        s2 = sinh(z/2)
        c2 = cosh(z/2)

        U1 = 2s2*c2/sqrt(-b)
        U2 = -2s2*s2/b
        U0 = 1.0 - b*U2
        U3 = (s - U1)/b
        U4 = (s^2/2.0 - U2)/b
        U5 = (s^3/6.0 - U3)/b

        return U0, U1, U2, U3, U4, U5
    else
        return 1.0, s, (s^2)/2, (s^3)/6, (s^4)/24, (s^5)/120
    end
end

# function universal_series(s, b)

# end