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

# tests from Vallado
T = Float64
pos = T.([1131.340, -2282.343, 6672.423])
vel = T.([-5.64305, 4.30333, 2.42879])
dt  = T(40.0*60.0)
gm  = T(398600.4415)# per vallado
# gm = 398600.0

chi = norm(vel)^2/2 - gm/norm(pos)
a   = -gm/(2.0*chi)
alp = 1/a
alp = 2/sqrt(dot(pos, pos)) - dot(vel, vel)/gm
# a   = 1/alp

posf, velf = Kepler.propagate(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)

# comapre with spice
vcat(posf, velf) - statef

# angular momentum conservation
cross(pos, vel)
cross(posf, velf)

cross(pos, vel) - cross(posf, velf)

cross(statef[1:3], statef[4:6]) - cross(pos, vel)

#mechanical energy conservation
((norm(velf)^2)/2 - gm/norm(posf)) - ((norm(vel)^2)/2 - gm/norm(pos))
((norm(statef[4:6])^2)/2 - gm/norm(statef[1:3])) - ((norm(vel)^2)/2 - gm/norm(pos))

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

posf, velf = Kepler.solve(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)

# comapre with spice
vcat(posf, velf) - statef

# angular momentum conservation
cross(pos, vel)
cross(posf, velf)

cross(pos, vel) - cross(posf, velf)

cross(statef[1:3], statef[4:6]) - cross(pos, vel)

#mechanical energy conservation
((norm(velf)^2)/2 - gm/norm(posf)) - ((norm(vel)^2)/2 - gm/norm(pos))
((norm(statef[4:6])^2)/2 - gm/norm(statef[1:3])) - ((norm(vel)^2)/2 - gm/norm(pos))

# try non-dimensionalized
DU   = norm(pos)
TU   = sqrt(DU^3/gm)
pos0 = pos/DU
vel0 = vel/DU*TU
dt0  = dt/TU
posf, velf = Kepler.solve(pos0, vel0, dt0, 1.0)
vcat(posf, velf) - SPICE.prop2b(1.0, [pos0..., vel0...], dt0)

# s = 0.4957408029292113
# s = 0.0

# # test a parabolic orbit
pos = [1.0, 0.0, 0.0]
vel = [0.0, (1.0 + 0e-11)*sqrt(2), 0.0]
# dt  =  1.234750019852072
dt = 0.001
gm = 1.0 

posf, velf = Kepler.solve(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)
vcat(posf, velf) - statef

cross(pos, vel)
cross(posf, velf)

cross(pos, vel) - cross(posf, velf)

cross(statef[1:3], statef[4:6]) - cross(pos, vel)

((norm(velf)^2)/2 - gm/norm(posf)) - ((norm(vel)^2)/2 - gm/norm(pos))
((norm(statef[4:6])^2)/2 - gm/norm(statef[1:3])) - ((norm(vel)^2)/2 - gm/norm(pos))

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