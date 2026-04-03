using Pkg; Pkg.activate()
using Kepler
using ForwardDiff
using StaticArrays

gm = 1.0
dt = pi/4

pos = SVector{3}(1.0, 0.0, 0.0)
vel = SVector{3}(0.0, 1.2, 0.0)

posf, velf, stm = Kepler.propagate_stm(pos, vel, dt, gm)

# check partials
dxdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[1], pos)
stm.dX_dX0 .- dxdx_auto

dxdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[1], vel)
stm.dX_dV0 .- dxdv_auto

dvdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[2], pos)
stm.dV_dX0 .- dvdx_auto

dvdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[2], vel)
stm.dV_dV0 .- dvdv_auto