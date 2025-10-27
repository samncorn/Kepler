# %%
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
using Roots
# %%

debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger)

# test from vallado
# T = Float64
# pos = T.([1131.340, -2282.343, 6672.423])
# vel = T.([-5.64305, 4.30333, 2.42879])
# dt  = T(40.0*60.0)
# gm  = T(398600.4415)# per vallado

# case that fails in the wild [FIXED]
# pos = [1.25, 0.0, 0.0]
# vel = [-0.021759125285088412, -0.015386025041773686, -0.0]
# dt  = 4.0
# gm  = 0.0002959122082326087

# parabolic orbit that failed with elliptic upper bound [FIXED]
pos = [1.25, 0.0, 0.0]
vel = [-0.010879562642544206, -0.018843955261014882, -0.0]
dt  = 8.0
gm  = 0.0002959122082326087

# # near rectilinear, consider as invalid input
# pos = [2.5, 0.0, 0.0]
# vel = [0.007693012520886843, 1e-11, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# failed in the wild
# runs fine, probably forgot to update to latest version 
# pos = [0.4186942098016838, 0.24637616481752644, 1.3672151662793137]
# vel = [0.023317918469691952, -0.014641183526584653, 0.0008003774012850409]
# dt  = -5.998977942803302
# gm  = 0.0002959122082326087

# test and compare to spice
posf, velf = Kepler.propagate(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)
posf2 = statef[1:3]
velf2 = statef[4:6]

posf - posf2
velf - velf2

# # break down internals
function universal_kepler(y, l, k1, k2, k3)
    _, c1, c2, c3 = Kepler.stumpff(l*y^2)
    L = y*(k1*c1 + y*(k2*c2 + y*k3*c3))
    return L
end

hvec = cross(pos, vel)
h2   = dot(hvec, hvec)
r0   = norm(pos)
dr0  = dot(vel, pos)/r0
v20  = dot(vel, vel)
a    = 2gm/r0 - v20
e    = sqrt(1 - a*h2/gm^2)
q    = h2/(gm*(1 + e))
l    = (1-e)/(1+e)

nu = sqrt(gm*(1 + e)/(q^3))
L  = nu*dt
k1 = r0/q
k2 = r0*dr0/(nu*(q^2))
k3 = gm/((nu^2)*(q^3))

bracket = (0.0, L)
x0 = sqrt(q/(1+e))*L
z0 = l*L^2 

universal_kepler(L, l, k1, k2, k3)/nu - dt

y = find_zero(_y -> universal_kepler(_y, l, k1, k2, k3) - L, bracket, A42())

# 
x = sqrt(q/(1+e))*y
l*y^2