using Kepler
using ForwardDiff
using StaticArrays
using SPICE
using LinearAlgebra

pos = SA[1131.340, -2282.343, 6672.423]
vel = SA[-5.64305, 4.30333, 2.42879]
dt  = 40.0*60.0
gm  = 398600.4415 # per vallado

# check partials
posf2, velf2, stm2 = Kepler.propagate_with_partials_old(pos, vel, dt, gm)
posf, velf, stm =  Kepler.propagate_with_partials(pos, vel, dt, gm)

posf - posf2
velf - velf2

stm.dX_dX0 - stm2.dX_dX0
stm.dX_dX0 ./ stm2.dX_dX0


dxdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[1], pos)
stm.dX_dX0 .- dxdx_auto
stm.dX_dX0 ./ dxdx_auto
stm.dX_dX0
dxdx_auto



# maximum(abs.(stm.dX_dX0 .- dxdx_auto)) 
# maximum(abs.(stm.dX_dX0 .- dxdx_auto) ./ stm.dX_dX0) 

dxdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[1], vel)
stm.dX_dV0 .- dxdv_auto 

# maximum(abs.(dxdv .- dxdv_auto))
# maximum(abs.(dxdv .- dxdv_auto) ./ dxdv)

dvdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[2], pos)
stm.dV_dX0 .- dvdx_auto 
# maximum(abs.(dvdx .- dvdx_auto) ./ dvdx)

dvdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[2], vel)
stm.dV_dV0 .- dvdv_auto