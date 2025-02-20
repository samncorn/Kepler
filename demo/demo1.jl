using UniversalVariable
using StaticArrays
using LinearAlgebra
using Printf
using SPICE

pos0 = @SVector [1.0, 0.0, 0.0]
vel0 = @SVector [0.0, 1.0, 0.0]
posf, velf = UniversalVariable.kepler_uv(pos0, vel0, pi/2, 1.0; max_iter = 10000)

for _ in 1:10
    pos = @SVector rand(3)
    vel = 2*(@SVector rand(3))
    t = rand()

    h = norm(cross(pos, vel))

    gm = 1.0

    posf, velf = UniversalVariable.kepler_uv(pos, vel, t, gm; anomaly_tol = 1e-8, max_iter = 1000)
    state = prop2b(gm, [pos..., vel...], t)
    posf2 = state[1:3]
    velf2 = state[4:6]

    hf  = norm(cross(posf, velf))
    hf2 = norm(cross(posf2, velf2))
    @printf "h error: %.3e\n" (hf - h)/h
    @printf "spice h error: %.3e\n" (hf2 - h)/h

    r0 = norm(pos)
    v02 = dot(vel, vel)
    b = 2gm/r0 - v02
    @printf "b: %.3f\n" b
    @printf "a: %.3f\n" gm/b

    @printf "pos: %.3e, %.3e, %.3e\n" (posf - posf2)...
    @printf "vel: %.3e, %.3e, %.3e\n" (velf - velf2)...
    println()
end
