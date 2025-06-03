using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using BenchmarkTools

pos1 = SA[1.0, 0.0, 0.0]
pos2 = SA[0.0, 1.0, 0.0]
dt = pi/2 
gm = 1.0

vel1, vel2 = Kepler.lambert_solve(pos1, pos2, dt, gm; max_iter = 10, tol = 1e-15)
posf, velf = Kepler.propagate(pos1, vel1, dt, gm)
posf - pos2
velf - vel2

# test some random orbits

pos1 = SA[1.0, 0.0, 0.0]
pos2 = SA[0.0, 1.5, 0.0]
dt = pi/2 
gm = 1.0

vel1, vel2 = Kepler.lambert_solve(pos1, pos2, dt, gm; max_iter = 100, tol = 1e-15)
posf, velf = Kepler.propagate(pos1, vel1, dt, gm)
posf - pos2
velf - vel2

@btime Kepler.lambert_solve($pos1, $pos2, $dt, $gm; max_iter = 100, tol = 1e-15)