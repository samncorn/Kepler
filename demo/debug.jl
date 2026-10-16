using LinearAlgebra
using StaticArrays
using Kepler
using BenchmarkTools
using SPICE
using Printf

function propagate(pos, vel, dt, gm)
    delp, delv = kepler_step(pos, vel, dt, gm)
    return pos + delp, vel + delv
end

function kepler_step(pos, vel, dt, gm)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)

    println()
    println("START")
    println()

    println("  r0   = $r0")
    println("  rdr0 = $rdr0")
    println("  b    = $b")

    # println("  e = ")

    s = solve_universal_kepler(dt, r0, rdr0, b, gm)
    f, g, df, dg = lagrange_step_coefficients(pos, vel, dt, gm, s)
            
    delp =  f*pos +  g*vel
    delv = df*pos + dg*vel
    return delp, delv
end

function lagrange_step_coefficients(pos, vel, dt, gm, s)
    r0   = norm(pos)
    rdr0 = dot(vel, pos)
    b    = 2gm/r0 - dot(vel, vel)
    
    U = Kepler.universal03(b, s) # index i -> i + 1
    f = -gm/r0*U[3]
    g = r0*U[2] + rdr0*U[3]
    
    posf = pos + f*pos + g*vel
    rf   = norm(posf)

    df = -gm/(r0*rf)*U[2]
    dg = -gm/rf*U[3]
    return f, g, df, dg
end

function solve_universal_kepler(dt, r0, rdr0, b, gm)
    # establish initial guess using a newton step from zero (what about a halley step?)
    sl = zero(dt)
    sh = dt/r0

    # make sure it's a bracket
    fl = -dt
    fh = universal_kepler(r0, rdr0, b, gm, sh) - dt
    println("  sl = $sl, fl = $fl")
    println("  sh = $sh, fh = $fh")
    i = 0
    while sign(fl)*sign(fh) > 0 || (isfinite(fl) && !isfinite(fh))
        i += 1
        if isfinite(fh)
            # undershot, shift the bracket
            println(  "...shifting")
            (sl, fl) = (sh, fh)
            sh *= 2
        else
            # overflow but lower bound still valid, bisect
            println(  "...squeezing")
            sh = sl + (sh - sl)/2
        end
        fh = universal_kepler(r0, rdr0, b, gm, sh) - dt
    end

    println("  bracket shifted $i times")
    println("  sl = $sl, fl = $fl")
    println("  sh = $sh, fh = $fh")

    if fl == 0
        return sl
    end

    if fh == 0
        return sh
    end

    @assert isfinite(fl) "Failed to find a finite bracket given input type $(typeof(fh))"

    # return chandrupatla(s -> universal_kepler(r0, rdr0, b, gm, s) - dt, (sl, sh))
    function univ_kep(s)
        t, r = universal_kepler_with_deriv(r0, rdr0, b, gm, s)
        return t - dt, r
    end
    return newton_bracket(univ_kep, (sl, sh))
end

function universal_kepler(r0, rdr0, b, gm, s)
    U = Kepler.universal03(b, s)
    t = r0*U[2] + rdr0*U[3] + gm*U[4]
    return t
end

function universal_kepler_with_deriv(r0, rdr0, b, gm, s)
    U = Kepler.universal03(b, s)
    t = r0*U[2] + rdr0*U[3] + gm*U[4]
    r = r0*U[1] + rdr0*U[2] + gm*U[3]
    return t, r
end

function chandrupatla(f::F, bracket) where {F}
    t = 0.5
    x2, x3 = bracket

    y2 = f(x2)
    y3 = f(x3)

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

        println()
        println(i)
        # println((x1 = x1, x2 = x2, x3 = x3))
        # println((y1 = y1, y2 = y2, y3 = y3))
        @printf "  x1 = %13.6e, x2 = %13.6e, x3 = %13.6e\n" x1 x2 x3
        @printf "  y1 = %13.6e, y2 = %13.6e, y3 = %13.6e\n" y1 y2 y3
        

        phi = (y1 - y2)/(y3 - y2)
        xi  = (x1 - x2)/(x3 - x2)

        if 4dx1 < dx3 && phi^2 < xi < 1 - (1 - phi)^2
            println("interpolate")
            a = (x3 - x1)/(x2 - x1)
            t = (y1/(y1 - y2))*(y3/(y3 - y2)) - a*(y1/(y3 - y1))*(y2/(y2 - y3))
            # t *= 0.95

            # if abs(y1) < 1e-13
            #     # safeguard against too good of a guess
            #     # t = min(0.5, 4*t)
            #     t *= 2
            # end
            
        else 
            println("bisect")
            t = 0.5
        end

        if t >= 1
            t = 0.5
        end
        
        # if abs(y1) < 1e-13
        #     # safeguard against too good of a guess
        #     # t = min(0.5, 4*t)
        #     t = 0.1
        # end
        
        (x3, y3, dx3) = (x2, y2, dx2)
        (x2, y2, dx2) = (x1, y1, dx1)
        x1  = x2 + t*(x3 - x2)
        dx1 = abs(t*(x3 - x2))
        # println((x2 = x2, x1 = x2 + t*(x3 - x2)))
    end
    println()
    println("  i = $i")
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

        println()
        println(i)
        @printf "  x = %13.6e, x2 = %13.6e, x1 = %13.6e\n" x x2 x1
        @printf "  y = %13.6e, y2 = %13.6e, y1 = %13.6e\n" y y2 y1

        _x = if 4dx1 < dx3
            println("newton")
            (x*y[2] - y[1])/y[2]
        else
            println(bisect)
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
    println()
    println("  i = $i")
    return x
end

# pos = [-1.0303099831796936, 0.14602260557381744, 0.0]
# vel = [-0.0023898026898215245, -0.016521419336222164, 0.0]
# dt = -1.826507356738252
# gm = 0.00029584
# 
gm = 0.0172^2

ei  = 0.0001
dti = 1e-6

pos = [1.0, 0.0, 0.0]
rp  = norm(pos)
a   = (1 - ei)/rp
# @info "  semimajor axis = $(a)"
# @info "eccentricity = $ei"
# p = 
# n   = 2sqrt(abs(p^(-3)))
n   = sqrt(gm)
T   = 2pi/n
v   = sqrt(2gm/rp - gm*a)
vel = [0.0, v, 0.0]
dt  = dti*T


hvec = cross(pos, vel)
evec = cross(vel, hvec)/gm - normalize(pos)
e = norm(evec)

a = 2/norm(pos) - dot(vel, vel)/gm
T = 2pi*sqrt(a^3/gm)
dt/T

posf, velf = propagate(pos, vel, dt, gm)
xs = SPICE.prop2b(gm, vcat(pos, vel), dt)

poss = xs[1:3]
vels = xs[4:6]

posf .- poss
velf .- vels

p = SA[pos...]
v = SA[vel...]
@btime Kepler.propagate($p, $v, $dt, $gm)
