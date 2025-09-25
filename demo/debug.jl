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

# case that fails in the wild [FIXED]
# pos = [1.25, 0.0, 0.0]
# vel = [-0.021759125285088412, -0.015386025041773686, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# parabolic orbit that failed with elliptic upper bound [FIXED]
# pos = [1.25, 0.0, 0.0]
# vel = [-0.010879562642544206, -0.018843955261014882, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# another parabolic failure
pos = [2.5, 0.0, 0.0]
vel = [0.007693012520886843, 1e-11, -0.0]
dt  = 8.0
gm  = 0.0002959122082326087

posf, velf = Kepler.propagate(pos, vel, dt, gm)

DU   = norm(pos)
TU   = sqrt(DU^3 / gm)
pos0 = pos/DU
vel0 = vel/DU*TU
dt0  = dt/TU

norm(vel0)^2

dr0 = dot(vel0, pos0)
a   = 2.0 - dot(vel0, vel0)

Hvec = cross(pos0, vel0)
    # Evec = cross(vel0, Hvec)/gm - pos0/r0
Evec = cross(vel0, Hvec) - pos0
e    = norm(Evec)

x = dt0*a/(1 - e)