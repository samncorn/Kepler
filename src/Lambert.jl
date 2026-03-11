@enum LambertMotion Long Short
@enum LambertEnergy High Low

# """ Returns all possible solutions

# implementation of Izzo 2015
# """
# function lambert_solve(pos1, pos2, dt, gm;)
#     ,
# end

""" Returns only a single solution, if it exists

Otherwise returns nothing

Modification of the implementation from Izzo 2015

Note that we follow Izzo's paper in assuming the short direction of travel
"""
function lambert_solve(pos1, pos2, dt, gm, n_revs; long_arc = false, high_energy = false)
    @assert t > 0
    @assert g > 0
    c_vec = pos2 - pos1
    c  = norm(c_vec)
    r1 = norm(pos1)
    r2 = norm(pos2)
    s  = (r1 + r2 + c)/2
    r1_hat = pos1 ./ r1
    r2_hat = pos2 ./ r2
    h_hat  = normalize(cross(r1hat, r2hat))
    l2 = 1 - c/s
    l  = sqrt(l2)
    if long_arc
        l *= -1
    end

    t1_hat, t2_hat = if pos1[1]*pos2[2] - pos1[2]*pos2[1] < 0
        l *= -1
        (
            cross(r1_hat, h_hat),
            cross(r2_hat, h_hat)
        )
    else
        (
            cross(r1_hat, h_hat),
            cross(r2_hat, h_hat)
        )
    end

    T  = sqrt(2*gm/s^3)*dt
    xy = lambert_xy(l, T, n_revs; high_energy = high_energy)
    if isnothing(xy)
        # no valid solution
        return nothing
    end
    (x, y) = xy

    # reconstruct velocity vectors
    g = sqrt(gm*s/2)
    p = (r1 - r2)/c
    o = sqrt(1 - p^2)
    v1r =  g*((l*y - x) - p*(l*y + x))/r1
    v2r = -g*((l*y - x) - p*(l*y + x))/r2
    v1t =  g*o*(y + l*x)/r1
    v2t =  g*o*(y + l*x)/r2
    vel1 = v1r*r1_hat + v1t*t1_hat
    vel2 = v2r*r2_hat + v2t*t2_hat 
    return vel1, vel2
end

function lambert_xy(l, T, n_revs; high_energy = false)
    # iterate
    T0 = acos(l) + l*sqrt(1 - l^2)
    # get initial guess
    x0 = if n_revs == 0
        # single solution
        T1 = 2*(1 - l^3)/3
        x0 = if T >= T0
            (T0/T)^(2/3) - 1.0
        elseif T < T1
            (5/2)*T1*(T1 - T)/T/(1 - l^5) + 1.0
        else
            (T0/T)^log2(T1/T0) - 1.0
        end
        x0
    else
        # find T_min
        T0 += n_revs*pi
        if T < T00
            # check for  existence of solution
            T_min = find_zero(x -> lambert_tof_with_derivative(x, l, n_revs)[2], 0.0)
            if T < T_min
                return nothing
            end
        end

        k = if high_energy
            (((n_revs+1)*pi)/(8T))^(2/3)
        else
            ((8T)/(n_revs*pi))^(2/3)
        end
        x0 = (k-1)/(k+1)
        x0
    end

    # find value
    x = find_zero(x -> lambert_tof(x, l, n_revs) - T, x0)
    y = sqrt(1 - l^2*(1 - x^2))
    return (x, y)
end

# TODO: 
# returns all valid cases
# function lambert_xy(l, T)

# end

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

