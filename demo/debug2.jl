# %%
using Pkg; Pkg.activate()
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
using Roots
using BenchmarkTools
# using ForwardDiff
using CairoMakie

import Test: @inferred
using Profile
# %%

# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debug_logger)

# test from vallado
# pos0 = SA[1131.340, -2282.343, 6672.423]
# vel0 = SA[-5.64305, 4.30333, 2.42879]
# dt0  = 40.0*60.0
# gm   = 398600.4415 # per vallado
# pos0 = [0.9230171944187261, -0.39465472915751687, 0.0]
# vel0 = [0.0065990688167711494, 0.016273154705785362, 0.0]
# dt0 = 3.653014713476504
# gm = 0.00029584

pos0 = [-0.8823881190038564, 2.750761917003355, 0.0]
vel0 = [-0.011570108698449117, 0.008474907261080352, 0.0]
dt0  = -1.826507356738252
gm   = 0.00029584

DU = norm(pos0)
TU = sqrt(DU^3/abs(gm))

pos = pos0/DU
vel = vel0/DU*TU
dt  = dt0/TU

s0 = dot(vel, pos)
b  = 2.0 - dot(vel, vel)

x1 = zero(dt)
y1 = -dt
r1 = one(dt)
p1 = (x = x1, y = y1, dy = r1)

# initial guess
# x2     = Kepler.kepler_guess_canonical(pos, vel, dt)
x2     = dt
y2, r2 = Kepler.universal_kepler2_canonical(x2, b, s0)
y2    -= dt
p2 = (x = x2, y = y2, dy = r2)

# bracket
i = 0
while sign(p2.y) == sign(p1.y) && i < 1000
    i += 1
    if abs(p2.y) == Inf
        println("halving")
        # x2 /= 2
        x2 = p1.x + sign(dt)*(p2.x - p1.x)/2
        y2, r2 = Kepler.universal_kepler2_canonical(x2, b, s0)
        y2    -= dt
        p2 = (x = x2, y = y2, dy = r2)
    else
        println("stepping")
        # use a newton step to try to bracket
        p1 = p2

        # x2 = (x2*r2 - 2y2)/r2
        x2 *= 2
        y2, r2 = Kepler.universal_kepler2_canonical(x2, b, s0)
        y2    -= dt
        p2 = (x = x2, y = y2, dy = r2)
    end
end

p1
p2

x = Kepler.flmsm1_step(p1, p2)
if x < p1.x || x > p2.x
    println("bisecting")
    x = (p1.x + p2.x)/2
end
y, r = Kepler.universal_kepler2_canonical(x, b, s0)
y   -= dt
p3   = (x = x, y = y, dy = r)

if sign(p3.y) == sign(p2.y)
    println("shuffling")
    (p1, p2) = (p2, p1)
end

# i = 0
x = NaN # just it initialize the loop

v_x = []
for i in 1:10
    println(i)
    signab = sign((p1.y - p2.y)/(p1.x - p2.x))

    if signab == sign(p1.dy) && signab == sign(p2.dy) && signab == sign(p3.dy) && p1.y != p3.y
        println("interpolating")
        # derivatives are well suited for inverse interpolation
        if p2.y != p3.y
            # need unique values
            x = Kepler.flmsm1_step(p1, p2, p3)
        else
            x = Kepler.flmsm1_step(p1, p3)
        end
    else
        # attempt newton's method
        println("newton step")
        x = (p3.x*p3.dy - p3.y)/p3.dy
    end

    # if we have stepped outside the bracket, reflect
    if (x < p3.x && x < p2.x) || (x > p3.x && x > p2.x) # || abs((x - p3.x)/p3.x) < 2eps(p3.x)
        println("bisecting")
        # x = (p2.x + p3.x)/2
        x = 2(p3.x) - x
    end
    dx = x - p3.x
    println("p1 = $p1")
    println("p2 = $p2")
    println("p3 = $p3")
    println("abs step = $dx")
    println("rel step = $(dx/p3.x)")
    println("   new x = $x")

    if isnan(x)
        break
    end

    push!(v_x, x)

    y, r = Kepler.universal_kepler2_canonical(x, b, s0)
    y   -= dt
    p    = (x = x, y = y, dy = r)
    println("   new y = $y")

    if x == p3.x || x == p2.x 
        break
    end

    p1 = p2
    p2 = p3
    p3 = p
    if sign(p3.y) == sign(p2.y)
        (p1, p2) = (p2, p1)
    end
end

# Kepler.universal_kepler_canonical(prevfloat(prevfloat(x)), b, s0) - dt

f = Figure()
ax = Axis(f[1, 1])

t = 0:(x2/100):3
l = Kepler.universal_kepler_canonical.(t, b, s0)

lines!(ax, t, l)
for xi in v_x
    yi = Kepler.universal_kepler_canonical.(xi, b, s0)
    scatter!(ax, xi, yi)
end

w = 1e-14
xlims!(ax, v_x[end] - w, v_x[end] + w)

f