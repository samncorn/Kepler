using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using BenchmarkTools
using ProgressMeter

@showprogress for i in 1:100
    pos1 = SA[1.0, 0.0, 0.0]
    vel1 = SA[rand(3)...]
    # dt = rand()*2pi
    gm = 1.0

    a = 1/(2/norm(pos1) - dot(vel1, vel1)/gm)

    P = if a > 0
        2pi*sqrt(a^3/gm)
    else
        rand()
    end
    dt = rand()*P

    hvec = cross(pos1, vel1)
    evec = cross(vel1, hvec)/gm - normalize(pos1)
    e = norm(evec)

    # @printf "%d %.3e %.3e" i a e

    posf, velf = Kepler.propagate(pos1, vel1, dt, gm)

    vel1_l, velf_l = Kepler.lambert_direct(pos1, posf, dt, gm)

    # println(vel1_l - vel1)
    # println(velf_l - velf)
    @printf "dt = %.3e, a = %.3e, |dv1| = %.3e, |dv2| = %.3e\n" dt a norm(vel1_l - vel1) norm(velf_l - velf)
end


pos1 = SA[1.0, 0.0, 0.0]
vel1 = SA[0.0, 0.0, 0.0]

gm = 1.0

a = 1/(2/norm(pos1) - dot(vel1, vel1)/gm)

hvec = cross(pos1, vel1)
evec = cross(vel1, hvec)/gm - normalize(pos1)
e = norm(evec)

# unroll
begin
pos1 = SA[1.0, 0.0, 0.0]
vel1 = SA[rand(3)...]
# dt = rand()*2pi/10
gm = 1.0

a = 1/(2/norm(pos1) - dot(vel1, vel1)/gm)

hvec = cross(pos1, vel1)
evec = cross(vel1, hvec)/gm - normalize(pos1)
e = norm(evec)

P = 2pi*sqrt(a^3/gm)

dt = rand()*P

@printf "%.3e %.3e %.3e\n" a e dt/P

pos2, vel2 = Kepler.propagate(pos1, vel1, dt, gm)

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
# x, y = Kepler.lambert_xy_0rev(l, T)
# unroll the xy computation
T0 = acos(l) + l*sqrt(1 - l^2)
T1 = 2*(1 - l^3)/3
x  = if T >= T0
    (T0/T)^(2/3) - 1.0
elseif T < T1
    (5/2)*T1*(T1 - T)/T/(1 - l^5) + 1.0
else
    # ((T0/T)^log(T1/T0)) - 1.0
    exp(log(2) * log(T/T0) / log(T1/T0)) - 1
end

_T, _, _, _ = Kepler.lambert_tof_with_3_derivatives(x, l, 0)
T - _T

# xp = Inf
# dx = Inf
dx1 = NaN
dx  = NaN
i  = 0
while abs(dx1) != abs(dx) && i < 100
# for i in 1:100
    i += 1
    Ti, dT1, dT2, dT3 = Kepler.lambert_tof_with_3_derivatives(x, l, 0)
    Ti -= T
    dx  = Ti*(dT1^2 - Ti*dT2/2)/(dT1*(dT1^2 - Ti*dT2) + (dT3*Ti^2)/6) # householder 3 iterations
    xp  = x
    x  -= dx
    println((x, x + dx, Ti, dx))
end

_T, _, _, _ = Kepler.lambert_tof_with_3_derivatives(x, l, 0)
T - _T
y = sqrt(1 - l^2*(1 - x^2))
g = sqrt(gm*s/2)
p = (r1 - r2)/c
o = sqrt(1 - p^2)
v1r =  g*((l*y - x) - p*(l*y + x))/r1
v2r = -g*((l*y - x) + p*(l*y + x))/r2
v1t =  g*o*(y + l*x)/r1
v2t =  g*o*(y + l*x)/r2
vel1l = v1r*ir1 + v1t*it1
vel2l = v2r*ir2 + v2t*it2
end

norm(vel1 - vel1l)
norm(vel2 - vel2l)
@printf "%.3e %.3e %.3f\n" a e dt/P

using CairoMakie
f = Figure()
ax = Axis(f[1, 1])

a = -1:0.01:2
b = [Kepler.lambert_tof_with_3_derivatives(xi, l, 0)[1] for xi in a]
lines!(ax, a, b)
f