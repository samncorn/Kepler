function stumpff(z)
    if z > 0
        sin2 = sin(0.5sqrt(z))
        cos2 = cos(0.5sqrt(z))
        c1 = 2sin2*cos2/sqrt(z)
        c2 = 2sin2^2/z
        c3 = (1 - c1)/z
        c0 = 1 - z*c2
        return c0, c1, c2, c3
    elseif z < 0
        sinh2 = sinh(0.5sqrt(-z))
        cosh2 = cosh(0.5sqrt(-z))
        c1 = 2sinh2*cosh2/sqrt(-z)
        c2 = -2sinh2^2/z
        c3 = (1 - c1)/z
        c0 = 1 - z*c2
        return c0, c1, c2, c3
    else
        c1 = 1/2
        c2 = 1/6
        c3 = 1 - z*c1
        c0 = 1 - z*c2
        return c0, c1, c2, c3
    end
end
# function stumpff(b, s)
#     if b > 0
#         sin2 = sin(0.5sqrt(b)*s)
#         cos2 = cos(0.5sqrt(b)*s)

#         c1 = 2sin2*cos2/sqrt(b)
#         c2 = 2sin2*sin2/b
#         c3 = (s - c1)/b
#         c0 = 1 - b*c2
#         return c0, c1, c2, c3
#     elseif b > 0
#         sinh2 = sinh(0.5sqrt(-b)*s)
#         cosh2 = cosh(0.5sqrt(-b)*s)

#         c1 = 2sinh2*cosh2/sqrt(-b)
#         c2 = -2sinh2*sinh2/b
#         c3 = (s - c1)/b
#         c0 = 1 - b*c2
#         return c0, c1, c2, c3
#     else
#         c2 = 1/2
#         c3 = 1/6
#         c0 = 1.0
#         c1 = s
#         return c0, c1, c2, c3
#     end
# end

# function stumpff_series(b, s)

# end