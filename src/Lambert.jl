""" Returns only the direct (0-rev) solution.

Modification (barely) of the implementation from Izzo 2015
"""
function lambert_direct(pos1, pos2, dt, gm)
    c_vec = pos2 - pos1
    c  = norm(c_vec)
    r1 = norm(pos1)
    r2 = norm(pos2)
    s  = (r1 + r2 + c)/2
    ir1 = pos1 ./ r1
    ir2 = pos2 ./ r2
    ih  = normalize(cross(ir1, ir2))
    l2  = 1 - c/s
    l   = sqrt(l2)
    it1 = cross(ih, ir1)
    it2 = cross(ih, ir2)
    if pos1[1]*pos2[2] - pos1[2]*pos2[1] < 0
        l   *= -1
        it1 *= -1
        it2 *= -1
    end
    T    = sqrt(2gm/s^3)*dt
    x, y = lambert_xy_0rev(l, T)
    # reconstruct velocity vectors
    g = sqrt(gm*s/2)
    p = (r1 - r2)/c
    o = sqrt(1 - p^2)
    v1r =  g*((l*y - x) - p*(l*y + x))/r1
    v2r = -g*((l*y - x) + p*(l*y + x))/r2
    v1t =  g*o*(y + l*x)/r1
    v2t =  g*o*(y + l*x)/r2
    vel1 = v1r*ir1 + v1t*it1
    vel2 = v2r*ir2 + v2t*it2
    return vel1, vel2
end

function lambert_xy_0rev(l, T)
    T0 = acos(l) + l*sqrt(1 - l^2)
    T1 = 2*(1 - l^3)/3
    x  = if T >= T0
        (T0/T)^(2/3) - 1.0
    elseif T < T1
        (5/2)*T1*(T1 - T)/T/(1 - l^5) + 1.0
    else
        # (T0/T)^log2(T1/T0) - 1.0 # WRONG, PAPER IS WRONG
        exp(log(2) * log(T/T0) / log(T1/T0)) - 1
    end
    # iterate
    # need to bracket to avoid possible limit cycles
    # xp = Inf
    # dx = Inf
    dx1 = NaN
    dx  = NaN
    i  = 0
    while abs(dx1) != abs(dx) && i < 100
        i  += 1
        Ti, dT1, dT2, dT3 = lambert_tof_with_3_derivatives(x, l, 0)
        Ti -= T
        dx  = Ti*(dT1^2 - Ti*dT2/2)/(dT1*(dT1^2 - Ti*dT2) + (dT3*Ti^2)/6) # householder 3 iterations
        xp  = x
        x  -= dx
    end
    y = sqrt(1 - l^2*(1 - x^2))
    return x, y
end

function lambert_tof_with_3_derivatives(x, l, n)
    y = sqrt(1 - l^2*(1 - x^2))
    if n == 0 && abs(x - 1) < 1e-3
        # use hypergeometric series to evaluate
        e  = y - l*x
        S1 = (1.0 - l - x*e)/2
        Q  = 4*hypergeometric(S1)/3
        # NOTE THAT WE ARE TRUNCATING THE DERIVATIVE HERE (but we don't need it to be particularly precise)
        T  = e*(e^2 * Q + 4*l)/2
        dT1 = 2*(l^5 - 1)
        dT2 = 0.0 # assume derivative is 0 (approximate as linear derivs near 1 )
        dT3 = 0.0
        return (T, dT1, dT2, dT3)
    end

    phi = if x < 1
        c = x*y + l*(1 - x^2)
        s = (y - l*x)*sqrt(1 - x^2)
        atan(s, c)
    else
        c = x*y - l*(x^2 - 1)
        s = (y - l*x)*sqrt(x^2 - 1)
        atanh(s/c)
    end
    T   = ((phi + n*pi)/sqrt(abs(1 - x^2)) - x + l*y)/(1 - x^2)
    dT1 = (3T*x - 2.0 + 2.0*(l^3)*x/y)/(1 - x^2)
    dT2 = (3T + 5x*dT1 + 2*(1 - l^2)*(l^3)/(y^3))/(1 - x^2)
    dT3 = (7x*dT2 + 8dT1 - 6*(1 - l^2)*(l^5)*x/(y^5))/(1 - x^2)
    return (T, dT1, dT2, dT3)
end


# function lambert_guess(l, T, n_revs; high_energy = false)
#     T0 = acos(l) + l*sqrt(1 - l^2)
#     if n_revs == 0
#         # single solution
#         T1 = 2*(1 - l^3)/3
        
#         return x0
#     else
#         # find T_min
#         T0 += n_revs*pi
#         if T < T00
#             # check for  existence of solution
#             # T_min = find_zero(x -> lambert_tof_with_derivative(x, l, n_revs)[2], 0.0)
#             dx    = Inf
#             x_min = 0.0
#             while x_min - xd != x_min
#                 T_min, dT_min = lambert_tof_with_derivative(x_min, l, n_revs)
#                 dx = dT_min/T_min

#                 x_min -= dx
#             end
#             if T < T_min
#                 return nothing
#             end
#         end

#         k = if high_energy
#             (((n_revs+1)*pi)/(8T))^(2/3)
#         else
#             ((8T)/(n_revs*pi))^(2/3)
#         end
#         x0 = (k-1)/(k+1)
#         return x0
#     end
# end

function lambert_tof(x, l, M)
    y = sqrt(1 - l^2*(1 - x^2))
    if abs(x - z) < 1e-3
        # use hypergeometric series to evaluate
        e  = y - l*x
        S1 = (1.0 - l - x*e)/2
        Q  = 4*hypergeometric(S1)/3
        T  = e*(e^2 * Q + 4*l)/2
        return T
    end

    c = x*y + l*(1 - x^2)
    s = (y - l*x)*sqrt(1 - x^2)
    phi = if l >= 0
        atan(s, c)
    else
        atanh(s/c)
    end
    T  = ((phi + M*pi)/sqrt(abs(1 - x^2)) - x + l*y)/(1 - x^2)
    return T
end

function lambert_tof_with_derivative(x, l, M)
    y = sqrt(1 - l^2*(1 - x^2))
    if abs(x - 1) < 1e-3
        # use hypergeometric series to evaluate
        e  = y - l*x
        S1 = (1.0 - l - x*e)/2
        Q  = 4*hypergeometric(S1)/3
        # NOTE THAT WE ARE TRUNCATING THE DERIVATIVE HERE (but we don't need it to be particularly precise)
        T  = e*(e^2 * Q + 4*l)/2
        dt = 2*(l^5 - 1)
        return (T, dT)
    end

    c = x*y + l*(1 - x^2)
    s = (y - l*x)*sqrt(1 - x^2)
    phi = if l >= 0
        atan(s, c)
    else
        atanh(s/c)
    end
    T  = ((phi + M*pi)/sqrt(abs(1 - x^2)) - x + l*y)/(1 - x^2)
    dt = (3T*x - 2.0 + 2.0*(l^3)*x/y)/(1 - x^2)
    return (T, dT)
end

# hypergeometric series for a = 3, b = 1, c = 5/2
function hypergeometric(z)
    a  = 3.0
    b  = 1.0
    c  = 5.0/2.0
    x0 = 1.0
    xi = a*b*z/c
    x  = x0 + xi
    j  = 1.0
    while x != x0
        x0  = x
        xi *= (a + j)*(b + j)*z/(c + j)/(j + 1.0)
        x  += xi
    end
    return x
end

