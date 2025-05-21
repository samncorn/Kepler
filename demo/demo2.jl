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

# test a parabolic orbit
pos = [1.0, 0.0, 0.0]
vel = [0.0, (1.0 + 0e-11)*sqrt(2), 0.0]
# dt  =  1.234750019852072
dt = 0.001
gm = 1.0 

posf, velf = Kepler.solve(pos, vel, dt, gm; parabolic_tol = 1e-6)
vcat(posf, velf) - SPICE.prop2b(gm, [pos..., vel...], dt)
