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

    # if b == 0
    #     return 1.0, s, (s^2)/2, (s^3)/6
    # end

    # z = sqrt(b*s^2)
    # s2, c2 = if b < 0
    #     (
    #         sinh(z/2),
    #         cosh(z/2),
    #     )
    # else
    #     (
    #         sin(z/2),
    #         cos(z/2),
    #     )
    # end

    # U1 = 2s2*c2/sqrt(abs(b))
    # U2 = 2s2*s2/abs(b)
    # U0 = 1.0 - b*U2
    # U3 = (s - U1)/b

    # return U0, U1, U2, U3
    
end

function universal05(b, s)
    if b > 0
        z  = sqrt(b)*s
        s2 = sin(z/2)
        c2 = cos(z/2)
        # c2 = sqrt(1.0 - s2^2)

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
        # c2 = sqrt(1.0 + s2^2)

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


# function universal_battin(w, )
#     z = 
# end

# function solve_kepler_battin(dt, r0, s0, b, gm)
#     T = sqrt(gm/r0^3)*dt
#     # constants
#     psi0 = 
#     gam0 =

#     # initial guess
#     phi0 = T - psi0*T^2/2 - (1 - gam0 - 3psi0^2)*T^3/6 + psi0*(10 - 9*gam0 - 15*psi0^2)*T^4/24
#     w0   = if gam0 > 0
#         u0 = 
#         u1 =
#         u1/u0
#     elseif gam0 < 0
#         u0 = 
#         u1 =
#         u1/u0
#     else
#         phi0
#     end

    

# end