using Pkg; Pkg.activate()
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using ProgressMeter
using ForwardDiff
using GLMakie

# using au, au/day
const au = 149597870.7
const gm = 0.01720209894846^2

earth_kep = (
    a  = 1.4765067E+08/au,
    e  = 9.1669995E-03,
    i  = deg2rad(4.2422693E-03),
    w  = deg2rad(6.64375167E+01),
    Om = deg2rad(1.4760836E+01)
)

apophis_kep = (
    a  = 1.3793939E+08/au,
    e  = 1.9097084E-01,
    i  = deg2rad(3.3356539E+00),
    w  = deg2rad(1.2919949E+02),
    Om = deg2rad(2.0381969E+02)
)

yr4_kep = (
    a  = 3.7680703E+08/au,
    e  = 6.6164147E-01,
    i  = deg2rad(3.4001497E+00),
    w  = deg2rad(1.3429905E+02),
    Om = deg2rad(2.7147904E+02)
)

atlas_3I_kep = (
    a  = -3.9552667E+07/au,
    e  = 6.1469268E+00,
    i  = deg2rad(1.7512507E+02),
    w  = deg2rad(1.2817255E+02),
    Om = deg2rad(3.2228906E+02)
)

earth_com = (
    q  = earth_kep.a * (1.0 - earth_kep.e),
    e  = earth_kep.e,
    i  = earth_kep.i,
    w  = earth_kep.w,
    Om = earth_kep.Om,
)

apophis_com = (
    q  = apophis_kep.a * (1.0 - apophis_kep.e),
    e  = apophis_kep.e,
    i  = apophis_kep.i,
    w  = apophis_kep.w,
    Om = apophis_kep.Om,
)

# earth_com = (
#     q  = earth_kep.a * (1.0 - earth_kep.e),
#     e  = earth_kep.e,
#     i  = earth_kep.i,
#     w  = earth_kep.w,
#     Om = earth_kep.Om,
# )

# earth_com = (
#     q  = earth_kep.a * (1.0 - earth_kep.e),
#     e  = earth_kep.e,
#     i  = earth_kep.i,
#     w  = earth_kep.w,
#     Om = earth_kep.Om,
# )

earth_pos, earth_vel     = Kepler.cartesian(earth_com.q, earth_com.e, earth_com.i, earth_com.Om, earth_com.w, 0.0, 0.0, gm)
apophis_pos, apophis_vel = Kepler.cartesian(apophis_com.q, apophis_com.e, apophis_com.i, apophis_com.Om, apophis_com.w, 0.0, 0.0, gm)

pos1 = earth_pos
vel1 = earth_vel
p1   = norm(cross(pos1, vel1))^2/gm
e1   = earth_com.e
X1   = 

pos2 = apophis_pos
vel2 = apophis_vel
p2   = norm(cross(pos2, vel2))^2/gm
e2   = apophis_com.e

ellipse1 = Kepler.Ellipse()

# plot the cost space
# v_grid = 0:(2pi/100):2pi
v1_grid = (0:0.01:1) .* 2pi 
v2_grid = (0:0.01:1) .* 2pi 
M = zeros(Float64, length(v1_grid), length(v2_grid))
@showprogress for ((i1, v1), (i2, v2)) in Iterators.product(enumerate(v1_grid), enumerate(v2_grid))
    r1 = h12/gm/(1.0 + e1*cos(v1))
    x1 = r1*(cos(v1)*xhat1 + sin(v1)*yhat1)
    r2 = h22/gm/(1.0 + e2*cos(v2))
    x2 = r2*(cos(v2)*xhat2 + sin(v2)*yhat2)
    M[i1, i2] = norm(x1 - x2)
end

f  = Figure()
ax = Axis(f[1, 1])
hm = contourf!(ax, v1_grid, v2_grid, M; levels = 20)
f

# scanning pass
v2min = MVector{4}(Inf, Inf, Inf, Inf)
v1min = MVector{4}(Inf, Inf, Inf, Inf)
dmin = MVector{4}(Inf, Inf, Inf, Inf)
i    = 1
D1   = Inf
D2   = Inf
v1   = Inf
v2   = Inf
for v in 0:0.12:(2pi + 0.12)
    # find the position of o2
    r2   = h22/gm/(1.0 + e2*cos(v))
    pos2 = r2*(cos(v)*xhat2 + sin(v)*yhat2)
    
    # get out of plane componet w.r.t o1
    x2  = dot(pos2, xhat1)
    y2  = dot(pos2, yhat1)
    z22 = r2^2 - x2^2 - y2^2

    # compute o1 radial distance for same longitude
    rho = sqrt(x2^2 + y2^2)
    r1  = h12/gm/(1.0 + e1*x2/rho)

    # compute meriodional distance
    D = sqrt(z22 + (rho - r1)^2)

    # check for local min
    if D1 < D2 > D
        v1i = atan(y2, x2)
        if v1i < 0
            v1i = 2pi + v1i
        end
        v2min[i] = v
        v1min[i] = v1i
        dmin[i] = D2
        i += 1
    end

    D1 = D2
    D2 = D
    v1 = v2
    v2 = v
end

f  = Figure()
ax = Axis(f[1, 1])
hm = contourf!(ax, v_grid, v_grid, M; levels = 20)
for (v1, v2) in zip(v1min, v2min)
    if v2 == Inf
        continue
    end
    scatter!(ax, v1, v2)
end
f

# check for need to extend scan
n_dmin = 0
for d in dmin
    if d < Inf
        n_dmin += 1
    end
end
if true # n_dmin <= 1
    println("using 4 starting points...")
    # try four evenly spread starting points
    v2min = MVector{4}(0.0, pi/2, pi, 3pi/2)
    for (i, v) in enumerate(v2min)
        # find the position of o2
        r2   = h22^2/gm/(1.0 + e2*cos(v))
        pos2 = r2*(cos(v)*xhat2 + sin(v)*yhat2)
        
        # get out of plane componet w.r.t o1
        x2  = dot(pos2, xhat1)
        y2  = dot(pos2, yhat1)
        v1i = atan(y2, x2)
        if v1i < 0
            v1i = 2pi + v1i
        end
        v1min[i] = v1i
    end
end

f  = Figure()
ax = Axis(f[1, 1])
hm = contourf!(ax, v_grid, v_grid, M; levels = 100)
for (v1, v2) in zip(v1min, v2min)
    if v2 == Inf
        continue
    end
    scatter!(ax, v1, v2)
end
f

# function dist2(v1, v2)

# end

# function moid_cost_partials(v12, params)

# end

# dist2(v1, v2) = norm(
#       r1*(cos(v1)*xhat1 + sin(v1)*yhat1)
#     - r2*(cos(v2)*xhat2 + sin(v2)*yhat2)
# )

# dist2_partials(v1, v2) = 
function dists_with_partials(v1, h1, e1, X1, Y1, v2, h2, e2, X2, Y2, gm)
    r1 = h1/gm/(1.0 + e1*cos(v1))
    r2 = h2/gm/(1.0 + e2*cos(v2))

    dr1 = h1*e1*sin(v1)/gm/(1.0 + e1*cos(v1))^2
    dr2 = h2*e2*sin(v2)/gm/(1.0 + e2*cos(v2))^2

    pos1 = r1*(cos(v1)*X1 + sin(v1)*Y1)
    pos2 = r2*(cos(v2)*X2 + sin(v2)*Y2)
    dx   = pos1 - pos2
    H    = hcat(
         dr1*pos1/r1 + r1*(-sin(v1)*X1 + cos(v1)*Y1),
        -dr2*pos2/r2 - r2*(-sin(v2)*X2 + cos(v2)*Y2)
    )
    return dx, -H
end

# v1 = v1min[1]
# v2 = v2min[1]
# dx, H = dists_with_partials(v1, h12, e1, xhat1, yhat1, v2, h22, e2, xhat2, yhat2, gm)
# d   = dot(dx, dx)/2
# # Htx = transpose(H)*dx
# # HtH = transpose(H)*H
# Ht = transpose(H)
# (dv1, dv2) = pinv(transpose(Ht*H) + 100*I)*Ht*dx*10000
# arrows2d!(ax, (v1, v2), (dv1, dv2))

# v1 += dv1
# v2 += dv2
# scatter!(ax, v1, v2)

# dx, H = dists_with_partials(v1, h12, e1, xhat1, yhat1, v2, h22, e2, xhat2, yhat2, gm)
# Ht = transpose(H)
# (dv1, dv2) = pinv(Ht*H + 100*I)*Ht*dx*10000
# arrows2d!(ax, (v1, v2), (dv1, dv2))

# f

sink = [[], [], [], []]
moid = Inf
for (si, (v1, v2)) in enumerate(zip(v1min, v2min))
    if isinf(v1) || isinf(v2)
        continue
    end
# begin
    # v1 = v1min[1]
    # v2 = v2min[1]
    # v1 = rand()*2pi
    # v2 = rand()*2pi
    # try least squares stepping (levenberg marquardt)
    # initialize
    mu_f = 5.0
    mu   = 0.0

    dx, H = dists_with_partials(v1, h12, e1, xhat1, yhat1, v2, h22, e2, xhat2, yhat2, gm)
    d   = dot(dx, dx)/2
    Ht  = transpose(H)

    i  = 0
    d0 = 2d

    push!(sink[si], (
        i = i, 
        x = (v1, v2), 
        dx = (0.0, 0.0), 
        cost = d, 
        mu = mu,
    ))
    while abs(d - d0) > 1e-8 && i < 100
        i += 1
        # (dv1, dv2) = -pinv(Ht*H + mu*I)*Ht*dx
        (dv1, dv2) = pinv(Ht*H + mu*I)*Ht*dx
        # (dv1, dv2) = -Ht*dx

        dx, H = dists_with_partials(v1 + dv1, h12, e1, xhat1, yhat1, v2 + dv2, h22, e2, xhat2, yhat2, gm)
        di  = dot(dx, dx)/2
        Ht  = transpose(H)

        if di >= 1.1d
            mu *= mu_f
        else
            mu /= mu_f
        end

        v1 += dv1
        v2 += dv2
        v1  = mod(v1, 2pi)
        v2  = mod(v2, 2pi)
        d0  = d
        d   = di

        push!(sink[si], (
            i = i, 
            x = (v1, v2), 
            dx = (dv1, dv2), 
            cost = d, 
            mu = mu,
        ))
    end

    if d < moid
        moid = d
    end
end

f  = Figure()
ax = Axis(f[1, 1])
hm = contourf!(ax, v_grid, v_grid, M; levels = 20)
# hm = contourf!(ax, v_grid, v_grid, reverse(M; dims = 2); levels = 20)
markers = [
    :circle,
    :utriangle,
    :dtriangle,
    :rect
]
for (v1, v2, sinki, marker) in zip(v1min, v2min, sink, markers)
    if v2 == Inf
        continue
    end
    scatter!(ax, v1, v2)
    l = [x[2] for x in sinki]
    t = [x[1]/length(sink) for x in sinki]
    scatter!(ax, l; color = 1.0 .- t, marker = marker)
end
f

# armellin test cases
orbit1 = Kepler.Orbit(
    Kepler.cartesian(
        1.0,
        0.0,
        0.0,
        0.0,
        deg2rad(16.0),
        0.0, 
        0.0,
        gm
    )..., 
    0.0, 
    gm
)
orbit2 = Kepler.Orbit(
    Kepler.cartesian(
        0.48, 
        0.6, 
        deg2rad(60.0), 
        0.0, 
        deg2rad(176.0), 
        0.0, 
        0.0, 
        gm
    )..., 
    0.0, 
    gm
)

dmin, v1min, v2min = Kepler.moid_scan_meridional(orbit1.position, orbit1.velocity, orbit2.position, orbit2.velocity, gm)

# d, v1, v2 = Kepler.moid_scan(orbit1, orbit2)

# ellipse1.x
# # Kepler.vec_angle()
# 180 - atand(ellipse1.x[2], ellipse1.x[1])

# decompose
ellipse1 = Kepler.Ellipse(orbit1)
ellipse2 = Kepler.Ellipse(orbit2)

orbit1.position
ellipse1.x

hvec = cross(orbit1.position, orbit1.velocity)
evec = cross(orbit1.velocity, hvec)/gm - normalize(orbit1.position)

v1_grid = (0:0.01:1) .* 2pi 
v2_grid = (0:0.01:1) .* 2pi 
M = zeros(Float64, length(v1_grid), length(v2_grid))
@showprogress for ((i1, _v1), (i2, _v2)) in Iterators.product(enumerate(v1_grid), enumerate(v2_grid))
    dx, _ = Kepler.dists_with_partials(
        _v1, 
        ellipse1.p, 
        ellipse1.e, 
        ellipse1.x, 
        ellipse1.y, 
        _v2, 
        ellipse2.p, 
        ellipse2.e, 
        ellipse2.x, 
        ellipse2.y, 
    )
    M[i1, i2] = dot(dx, dx)
end

f  = Figure()
ax = Axis(f[1, 1])
hm = contourf!(ax, v1_grid, v2_grid, M; levels = 10)
scatter!(ax, v1, v2, d)
f
