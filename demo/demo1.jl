using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
using BenchmarkTools
# using CairoMakie
# using GLMakie

# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debug_logger)

pos0 = @SVector [1.0, 0.0, 0.0]
vel0 = @SVector [0.0, 1.1, 0.0]
posf, velf = Kepler.solve(pos0, vel0, pi/2, 1.0)
@btime Kepler.solve($pos0, $vel0, pi/2, 1.0)

for _ in 1:50
    q = rand()
    e = rand()

    gm = 1.0
    pos = q*SVector{3, Float64}([1.0, 0.0, 0.0])
    vel = sqrt(gm*(e/q + 1))*SVector{3, Float64}([0.0, 1.0, 0.0])
    t = rand()
    
    h = norm(cross(pos, vel))

    posf, velf = Kepler.solve(pos, vel, t, gm; max_iter = 1_000, tol = 1e-15)
    state = prop2b(gm, [pos..., vel...], t)
    posf2 = state[1:3]
    velf2 = state[4:6]

    @info "eccentricty: $(e)"

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

failed = []
succeeded = []
for _ in 1:1_000_000
    q = 3rand()
    e = 3rand()

    gm = 1.0
    t = rand()

    pos = q*SVector{3, Float64}([1.0, 0.0, 0.0])
    vel = sqrt(gm*(e/q + 1))*SVector{3, Float64}([0.0, 1.0, 0.0])

    try
        posf, velf = Kepler.solve(pos, vel, t, gm; max_iter = 1_000)
        
        if any(isnan.(posf))
            push!(failed, (q, e))
        else
            push!(succeeded, (q, e))
        end
    catch err
        @info "q: $(q)"
        @info "e: $(e)"
        # @info ""
        throw(err)
    end
end

q1 = [x[1] for x in succeeded]
e1 = [x[2] for x in succeeded]

q2 = [x[1] for x in failed]
e2 = [x[2] for x in failed]

#======================#

f = Figure(size = (700, 700))
ax = Axis(f[1, 1], xlabel = "periapse", ylabel = "eccentricity")

scatter!(ax, q1, e1, color = :gray, markersize = 2.0)
# scatter!(ax, q2, e2, color = :red, markersize = 5.0)

qline = 0:0.1:1.0
lines!(ax, qline, 1 .- qline)

f

#======================#

f = Figure(size = (700, 700))
ax = Axis(f[1, 1], xlabel = "semimajor axis", ylabel = "eccentricity")

scatter!(ax, q1 ./ (1 .- e1), e1, color = :gray, markersize = 2.0)
scatter!(ax, q2 ./ (1 .- e2), e2, color = :red, markersize = 5.0)

xlims!(ax, -5, 5)
# ylims!(ax, 0, 5)

f

#======================#

f = Figure(size = (700, 700))
ax = Axis(f[1, 1], xlabel = "alpha", ylabel = "eccentricity")

scatter!(ax, (1 .- e1) ./ q1, e1, color = :gray, markersize = 2.0)
scatter!(ax, (1 .- e2) ./ q2, e2, color = :red, markersize = 5.0)

xlims!(ax, -10, 10)
# ylims!(ax, 0, 5)

f

#======================#

failed = []
succeeded = []
for _ in 1:100_000
    gm = 1.0
    t = rand()

    pos = @SVector rand(3)
    vel = 10*(@SVector rand(3))

    r0 = norm(pos)
    v02 = dot(vel, vel)
    alpha = 2/r0 - v02 / gm
    a = 1/alpha

    H = cross(pos, vel)
    E = cross(vel, H)/gm - pos / r0
    e = norm(E)
    q = a*(1 - e)

    try
        posf, velf = Kepler.solve(pos, vel, t, gm; anomaly_tol = 1e-8, max_iter = 100)
        
        if any(isnan.(posf))
            push!(failed, (q, e))
        else
            push!(succeeded, (q, e))
        end
    catch e
        @info "a: $(a)"
        # @info ""
        throw(e)
    end
end

q1 = [x[1] for x in succeeded]
e1 = [x[2] for x in succeeded]

q2 = [x[1] for x in failed]
e2 = [x[2] for x in failed]

# g(x, y) = x^2 + y^2

# function test()
#     g1(x) = g(x, 1.0)
#     dg1   = ForwardDiff.derivative(g1, 1.0)

#     g2(x) = g(x, 2.0)
#     dg2   = ForwardDiff.derivative(g2, 1.0)
#     return g2(1.0), dg2
# end
