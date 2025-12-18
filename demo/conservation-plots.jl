using Kepler
using StaticArrays
using CairoMakie
using LinearAlgebra
using ProgressMeter
using SPICE

# work in normalized coords
# vary
    # timestep
    # eccentricity

# for each combination, the timestep is used to advance the orbit repeatedly through a wide number of phases to get a good sample
e_min = 0.0
e_max = 2.0
e_step = 0.02
e = e_min:e_step:e_max
length(e)

# fractions of the period
dt_step = 0.005
dt_max  = 0.5
dt = dt_step:dt_step:dt_max
length(dt)

# store the angular momentum and energy conservation
M_h = zeros(length(e), length(dt))
M_E = zeros(length(e), length(dt))

y = (sqrt(5) - 1)/2
gm = 0.0172^2

function try_prop(pos, vel, dt, gm)
    try
        return Kepler.propagate(pos, vel, dt, gm)
        # statef = SPICE.prop2b(gm, [pos..., vel...], dt)
        return statef[1:3], statef[4:6]
    catch _
        throw((pos = pos, vel = vel, dt = dt, gm = gm))
    end
end

items = collect(Iterators.product(enumerate(e), enumerate(dt)))
edt = Iterators.product(e, dt)

# @showprogress Threads.@threads for ((i, ei), (j, dti)) in items
@showprogress for ((i, ei), (j, dti)) in items
    # construct the periapse state
    pos = SVector{3}(1.0, 0.0, 0.0)
    rp  = norm(pos)
    a   = (1 - ei)/rp
    # p = 
    # n   = 2sqrt(abs(p^(-3)))
    n   = sqrt(gm)
    T   = 2pi/n
    v   = sqrt(2gm/rp - gm*a)
    vel = SVector{3}(0.0, v, 0.0)

    # jump from periapse to initial state
    t0 = 0
    while t0 < T/2
        pos, vel = try_prop(pos, vel, T*dti, gm)
        t0 += T*dti
    end
    pos, vel = try_prop(pos, vel, T*dti*y, gm)

    # compute initial conserved quantities
    h0 = norm(cross(pos, vel))
    E0 = dot(vel, vel)/2 - gm/norm(pos)

    # if E0 == 0
    #     throw("AHHHHHH")
    # end

    # jump around for a while
    for _ in 1:100
        while t0 > -T/2
            pos, vel = try_prop(pos, vel, -T*dti, gm)
            t0 -= T*dti
        end
        pos, vel = try_prop(pos, vel, T*dti*y, gm)
        while t0 < T/2
            pos, vel = try_prop(pos, vel, T*dti, gm)
            t0 += T*dti
        end
        pos, vel = try_prop(pos, vel, T*dti*y, gm)
    end

    hf = norm(cross(pos, vel))
    Ef = dot(vel, vel)/2 - gm/norm(pos)

    M_h[i, j] = 2(hf - h0)/(hf + h0)
    M_E[i, j] = 2(Ef - E0)/(Ef + E0)
    # M_h[i, j] = hf - h0
    # M_E[i, j] = Ef - E0
end

M_h_sign = sign.(M_h)
M_E_sign = sign.(M_E)

M_h = abs.(M_h)
M_E = abs.(M_E)

begin
f  = Figure(size = (800, 800))
ax = [
    Axis(f[1, 1][1, 1], title = "|h| error", xlabel = "eccentricity", ylabel = "dt/T"),
    Axis(f[2, 1][1, 1], title = "E error", xlabel = "eccentricity", ylabel = "dt/T"),
    Axis(f[1, 2][1, 1], title = "|h| error sign", xlabel = "eccentricity", ylabel = "dt/T"),
    Axis(f[2, 2][1, 1], title = "E error sign", xlabel = "eccentricity", ylabel = "dt/T"),
    ]
hm1 = heatmap!(ax[1], e, dt, M_h; colorrange = (1e-18, 1e-12), colorscale = log10)
hm2 = heatmap!(ax[2], e, dt, M_E; colorrange = (1e-18, 1e-12), colorscale = log10)
hm3 = heatmap!(ax[3], e, dt, M_h_sign)
hm4 = heatmap!(ax[4], e, dt, M_E_sign)

Colorbar(f[1, 1][1, 2], hm1)
Colorbar(f[2, 1][1, 2], hm2)
Colorbar(f[1, 2][1, 2], hm3)
Colorbar(f[2, 2][1, 2], hm4)

f
end

save(joinpath(@__DIR__, "../plots/error-spice.png"), f)
# save(joinpath(@__DIR__, "../plots/error.png"), f)