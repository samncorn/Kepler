# %%
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
# %%

debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger)

# a failure case encountered using heliolinc. hyperoblic, large |s|, sinh(sqrt(|z|)) overflows
# turns out it was a negative time issue (and the overflowed initial guess). chooseing a better guess and reversing time fixed it
# %%
# pos = [2.5, 0.0, 0.0]
# vel = [-0.0185, 4.8805672928708536e-5, 0.0]
pos = [2.5, 0.0, 0.0]
vel = [-0.015, 0.003378725763145275, 0.0]
dt = 1.234750019852072
# 
# 0.0002959122082326087
# dt  = -11.212502314709127
gm  = 0.01720209894846^2
# %%

posf, velf = Kepler.solve(pos, vel, dt, gm; parabolic_tol = 1e-6)
vcat(posf, velf) - SPICE.prop2b(gm, [pos..., vel...], dt)

# s = 0.4957408029292113
# s = 0.0

# # test a parabolic orbit
# pos = [1.0, 0.0, 0.0]
# vel = [0.0, (1.0 + 0e-11)*sqrt(2), 0.0]
# # dt  =  1.234750019852072
# dt = 0.001
# gm = 1.0 

# posf, velf = Kepler.solve(pos, vel, dt, gm; parabolic_tol = 1e-6)
# vcat(posf, velf) - SPICE.prop2b(gm, [pos..., vel...], dt)

using CairoMakie

# universal_kepler2(x, alpha, r0, dr0, gm)
# gm = 1.0

r0  = norm(pos)
dr0 = dot(vel, pos) / r0
v02 = dot(vel, vel)
alpha = 2gm/r0 - v02

sf = 0.4957408029292113
s0 = begin
    Hvec = cross(pos, vel)
    Evec = cross(vel, Hvec)/gm - pos/r0
    e = norm(Evec)
    dt*alpha/(gm*(1 - e))
end

dt1, r1 = Kepler.universal_kepler2(0.0, alpha, r0, dr0, gm)
dt2, r2 = Kepler.universal_kepler2(s0, alpha, r0, dr0, gm)

dtf = Kepler.universal_kepler(sf, alpha, r0, dr0, gm)

spl = Kepler.Cubic_Hermite_Spline(dt1, dt2, 0.0, s0, 1/r1, 1/r2)
s1  = Kepler.interpolate(dtf, spl)

dt3, r3 = Kepler.universal_kepler2(s1, alpha, r0, dr0, gm)
spl2    = Kepler.Cubic_Hermite_Spline(dt1, dt3, 0.0, s1, 1/r1, 1/r3)

s2 = Kepler.interpolate(dtf, spl2)
dt4, r4 = Kepler.universal_kepler2(s2, alpha, r0, dr0, gm)
dt4 - dtf

s = 0:0.01:20.0
times = Kepler.universal_kepler.(s, alpha, r0, dr0, gm)

f  = Figure(size = (900, 700))
ax = Axis(f[1, 1])

lines!(ax, s, times .- dtf)
lines!(ax, [Kepler.interpolate(t, spl) for t in times], times .- dtf)
lines!(ax, [Kepler.interpolate(t, spl2) for t in times], times .- dtf)

scatter!(ax, s0, Kepler.universal_kepler(s0, alpha, r0, dr0, gm) .- dtf)
scatter!(ax, sf, 0.0)
scatter!(ax, s1, Kepler.universal_kepler(s1, alpha, r0, dr0, gm) .- dtf)

f