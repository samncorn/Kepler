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

# case that fails in the wild
# pos = [1.25, 0.0, 0.0]
# vel = [-0.021759125285088412, -0.015386025041773686, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# parabolic orbit that failed with elliptic upper bound
pos = [1.25, 0.0, 0.0]
vel = [-0.010879562642544206, -0.018843955261014882, -0.0]
dt  = 8.0
gm  = 0.0002959122082326087

posf, velf = Kepler.propagate(pos, vel, dt, gm)

