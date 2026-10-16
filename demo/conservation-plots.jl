using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Kepler
using StaticArrays
using CairoMakie
using LinearAlgebra
using ProgressMeter
using SPICE
using Logging
using Dates

# for each combination, the timestep is used to advance the orbit repeatedly through a wide number of phases to get a good sample
debug_logger = ConsoleLogger(stderr, Logging.Warn)
global_logger(debug_logger)

n_steps = 100

# e_min = 0.0
# e_max = 2.0
# e_step = (e_max - e_min)/n_steps
# e = e_min:e_step:e_max
# length(e)
e = 1 .- 10 .^ (0:-0.2:-8)

# fractions of the period
# dt_step = 0.002
# dt_max  = 0.5
# dt_step = dt_max/n_steps
# dt = dt_step:dt_step:dt_max

dt = 10.0 .^ (-4:0.1:-0.3)
length(dt)

# store the angular momentum and energy conservation
M_h = zeros(length(e), length(dt))
M_E = zeros(length(e), length(dt))
M_t = zeros(length(e), length(dt))

y = (sqrt(5) - 1)/2
gm = 0.0172^2

function try_prop(pos, vel, dt, gm)
    Kepler.propagate(pos, vel, dt, gm)
#     try
#         # println((pos = pos, vel = vel, dt = dt, gm = gm))
#         # 
#         DU = 1.0
#         TU = 1.0
#         gmu = 1.0
# 
#         # DU = norm(pos)
#         # TU = sqrt(DU^3/abs(gm))
#         # gmu = gm
# 
#         posf, velf = Kepler.propagate(pos/DU, vel/DU*TU, dt/TU, gm/gmu)
#         
#         return posf*DU, velf*DU/TU
#         # statef = SPICE.prop2b(gm, [pos..., vel...], dt)
#         # return statef[1:3], statef[4:6]
#     catch _
#         throw((pos = pos, vel = vel, dt = dt, gm = gm))
#     end
end

items = collect(Iterators.product(enumerate(e), enumerate(dt)))
edt = Iterators.product(e, dt)

# @showprogress Threads.@threads for ((i, ei), (j, dti)) in item
@showprogress Threads.@threads for (ii, ((i, ei), (j, dti))) in collect(enumerate(items))
    t_start = time_ns()
    # println(ii)
    count = 1

    # @info "iteration $(ii)/$(length(items))"
    # construct the periapse state
    pos = SVector{3}(1.0, 0.0, 0.0)
    rp  = norm(pos)
    a   = (1 - ei)/rp
    # @info "  semimajor axis = $(a)"
    # @info "eccentricity = $ei"
    # p = 
    # n   = 2sqrt(abs(p^(-3)))
    n   = sqrt(gm)
    T   = 2pi/n
    v   = sqrt(2gm/rp - gm*a)
    vel = SVector{3}(0.0, v, 0.0)

    # jump from periapse to initial state
    t0 = 0.0
    jj = 0
    while t0 < T/2
        # println((ii, jj))
        # println((pos = pos, vel = vel, dt = T*dti, gm = gm))
        # @info "subiterate $(jj)"
        jj += 1
        pos, vel = try_prop(pos, vel, T*dti, gm)
        t0 += T*dti
        count += 1
        # 
    end
    # println((pos = pos, vel = vel, dt = T*dti*y, gm = gm))
    pos, vel = try_prop(pos, vel, T*dti*y, gm)
    count += 1

    # compute initial conserved quantities
    h0 = norm(cross(pos, vel))
    E0 = dot(vel, vel)/2 - gm/norm(pos)

    # if E0 == 0
    #     throw("AHHHHHH")
    # end

    # jump around for a while
    for _ in 1:100
        # @info "sweep 1"
        while t0 > -T/2
            # println((ii, jj))
            # println((pos = pos, vel = vel, dt = -T*dti, gm = gm))
            jj += 1
            # @info "  subiterate $(jj)"
            # @info "  pos = $pos"
            # @info "  vel = $vel"
            # @info "  dt  = $(-T*dti)"
            # @info "  gm  = $gm"
            pos, vel = try_prop(pos, vel, -T*dti, gm)
            count += 1
            t0 -= T*dti
            # 
        end
        pos, vel = try_prop(pos, vel, T*dti*y, gm)
        count += 1
        # println((pos = pos, vel = vel, dt = T*dti*y, gm = gm))
        while t0 < T/2
            # println((ii, jj))
            # println((pos = pos, vel = vel, dt = T*dti, gm = gm))
            jj += 1 
            # @info "subiterate $(jj)"
            pos, vel = try_prop(pos, vel, T*dti, gm)
            count += 1
            t0 += T*dti
        end
        # println((pos = pos, vel = vel, dt = T*dti, gm = gm))
        pos, vel = try_prop(pos, vel, T*dti*y, gm)
        count += 1
    end

    hf = norm(cross(pos, vel))
    Ef = dot(vel, vel)/2 - gm/norm(pos)

    # M_h[i, j] = (hf - h0)/h0 
    M_h[i, j] = (hf - h0)/(h0 + eps())
    # M_E[i, j] = abs(E0) > 1e-14 ? (Ef - E0)/E0 : Ef - E0
    M_E[i, j] = (Ef - E0)/E0
    
    # M_h[i, j] = hf - h0
    # M_E[i, j] = Ef - E0
    t_stop = time_ns()
    M_t[i, j] = (t_stop - t_start) / count / 1000
end

M_h_sign = sign.(M_h)
M_E_sign = sign.(M_E)

M_h = abs.(M_h)
M_E = abs.(M_E)

begin
f  = Figure(size = (800, 800))
ax = [
    Axis(f[1, 1][1, 1], title = "|h| error",        ylabel = "log_10(1 - e)", xlabel = "log10(dt/T)"),
    Axis(f[2, 1][1, 1], title = "|E| error",          ylabel = "log_10(1 - e)", xlabel = "log10(dt/T)"),
    Axis(f[1, 2][1, 1], title = "|h| error sign",   ylabel = "log_10(1 - e)", xlabel = "log10(dt/T)"),
    Axis(f[2, 2][1, 1], title = "E error sign",     ylabel = "log_10(1 - e)", xlabel = "log10(dt/T)"),
    # Axis(f[3, 1][1, 1], title = "run time",     ylabel = "log_10(1 - 2)", xlabel = "log10(dt/T)"),
    # Axis(f[3, 2][1, 1], title = "run time",     ylabel = "log_10(1 - 2)", xlabel = "log10(dt/T)"),
    ]

for axi in ax
    axi.yreversed = true
end

hm1 = heatmap!(ax[1], log10.(dt), log10.(1 .- e), transpose(M_h); colormap = :cubehelix, colorrange = (1e-16, 1e-6), colorscale = log10)
hm2 = heatmap!(ax[2], log10.(dt), log10.(1 .- e), transpose(M_E); colormap = :cubehelix, colorrange = (1e-16, 1e-6), colorscale = log10)
hm3 = heatmap!(ax[3], log10.(dt), log10.(1 .- e), transpose(M_h_sign))
hm4 = heatmap!(ax[4], log10.(dt), log10.(1 .- e), transpose(M_E_sign))
# hm5 = heatmap!(ax[5], log10.(dt), log10.(1 .- e), transpose(M_t))

Colorbar(f[1, 1][1, 2], hm1)
Colorbar(f[2, 1][1, 2], hm2)
Colorbar(f[1, 2][1, 2], hm3)
Colorbar(f[2, 2][1, 2], hm4)
# Colorbar(f[3, 1][1, 2], hm5)

f
end

# save(joinpath(@__DIR__, "../plots/error-spice.png"), f)
save(joinpath(@__DIR__, "../plots/error.png"), f)

# M_t[M_t .== 0] .= eps()

f  = Figure(size = (400, 400))
ax = Axis(f[1, 1], title = " average runtime (microsec)", ylabel = "log_10(1 - e)", xlabel = "log10(dt/T)")
hm = heatmap!(ax, log10.(dt), log10.(1 .- e), transpose(M_t); colorrange = (0, 2), colormap = :cubehelix)
Colorbar(f[1, 2], hm)
save(joinpath(@__DIR__, "../plots/runtime.png"), f)
