using Pkg; Pkg.activate()
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using ProgressMeter
using ForwardDiff
using CSV
using DataFrames
using Printf
using NonlinearSolve
using StatsBase
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

# given
apophis_kep = (
    a  = 1.3793939E+08/au,
    e  = 1.9097084E-01,
    i  = deg2rad(3.3356539E+00),
    w  = deg2rad(1.2919949E+02),
    Om = deg2rad(2.0381969E+02)
)

# jpls sbdb
# apophis_kep = (
#     a  = 0.9223803173917017,
#     e  = 0.1911663355386932,
#     i  = deg2rad(3.340958441017069),
#     w  = deg2rad(126.6728325163065),
#     Om = deg2rad(203.8996515621043)
# )

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

yr4_com = (
    q  = yr4_kep.a * (1.0 - yr4_kep.e),
    e  = yr4_kep.e,
    i  = yr4_kep.i,
    w  = yr4_kep.w,
    Om = yr4_kep.Om,
)

atlas_3I_com = (
    q  = atlas_3I_kep.a * (1.0 - atlas_3I_kep.e),
    e  = atlas_3I_kep.e,
    i  = atlas_3I_kep.i,
    w  = atlas_3I_kep.w,
    Om = atlas_3I_kep.Om,
)

apophis_com = (
    q  = apophis_kep.a * (1.0 - apophis_kep.e),
    e  = apophis_kep.e,
    i  = apophis_kep.i,
    w  = apophis_kep.w,
    Om = apophis_kep.Om,
)

apophis_eq = (
    a  =  0.922438242375914,
    P1 = −0.093144699837425,
    P2 =  0.166982492089134,
    Q1 =  0.012032857685451,
    Q2 = −0.026474053361345,
    l  =  deg2rad(88.3150906433494)
)

apophis_com2 = (
    e  = sqrt(apophis_eq.P1^2 + apophis_eq.P2^2),
    q  = (1 - sqrt(apophis_eq.P1^2 + apophis_eq.P2^2)) * apophis_eq.a,
    i  = mod(2atan(sqrt(apophis_eq.Q1^2 + apophis_eq.Q2^2)), 2pi),
    Om = mod(atan(apophis_eq.Q1/apophis_eq.Q2), 2pi),
    w  = mod(atan(apophis_eq.P1/apophis_eq.P2) - atan(apophis_eq.Q1/apophis_eq.Q2), 2pi),
)

cases = [
    # armellin test cases
    (tag = "Armellin 1", rp1 = 1.0,   e1 = 0.0,   rp2 = 0.48,  e2 = 0.6,   i1 = 0.0, i2 = 60.0, w1 = 16.0, w2 = 176.0, Om1 = 0.0, Om2 = 0.0, moid = 0.519407),
    (tag = "Armellin 2", rp1 = 0.585, e1 = 0.415, rp2 = 0.462, e2 = 0.615, i1 = 0.0, i2 = 80.0, w1 = 8.0,  w2 = 176.0, Om1 = 0.0, Om2 = 0.0, moid = 0.833579),
    (tag = "Armellin 3", rp1 = 1.0,   e1 = 0.6,   rp2 = 1.2,   e2 = 1.1,   i1 = 0.0, i2 = 40.0, w1 = 73.0, w2 = 69.0,  Om1 = 0.0, Om2 = 0.0, moid = 0.346197),
    (tag = "Armellin 4", rp1 = 1.0,   e1 = 0.5,   rp2 = 1.2,   e2 = 1.1,   i1 = 0.0, i2 = 66.0, w1 = 4.0,  w2 = 136.0, Om1 = 0.0, Om2 = 0.0, moid = 1.442149),
    # (rp1 = earth_com.q, e1 = earth_com.e, rp2 = apophis_com2.q, e2 = apophis_com2.e, i1 = rad2deg(earth_com.i), i2 = apophis_com2.i, w1 = rad2deg(earth_com.w), w2 = apophis_com2.w)
    # (rp1 = 1.0, e1 = 0.0, rp2 = apophis_com2.q, e2 = apophis_com2.e, i1 = 0.0, i2 = apophis_com2.i, w1 = 1.0, w2 = apophis_com2.w)
    
    # earth-apophis (hw)
    (
        tag = "Earth-Apophis",
        rp1 = earth_com.q, 
        e1  = earth_com.e, 
        rp2 = apophis_com.q, 
        e2  = apophis_com.e, 
        i1  = rad2deg(earth_com.i), 
        i2  = rad2deg(apophis_com.i), 
        w1  = rad2deg(earth_com.w), 
        w2  = rad2deg(apophis_com.w),
        Om1 = rad2deg(earth_com.Om), 
        Om2 = rad2deg(apophis_com.Om), 
        moid = NaN,
    ),
    (
        tag = "Earth-YR4",
        rp1 = earth_com.q, 
        e1  = earth_com.e, 
        rp2 = yr4_com.q, 
        e2  = yr4_com.e, 
        i1  = rad2deg(earth_com.i), 
        i2  = rad2deg(yr4_com.i), 
        w1  = rad2deg(earth_com.w), 
        w2  = rad2deg(yr4_com.w),
        Om1 = rad2deg(earth_com.Om), 
        Om2 = rad2deg(yr4_com.Om), 
        moid = NaN,
    ),
    (
        tag = "Apophis-YR4",
        rp1 = apophis_com.q, 
        e1  = apophis_com.e, 
        rp2 = yr4_com.q, 
        e2  = yr4_com.e, 
        i1  = rad2deg(apophis_com.i), 
        i2  = rad2deg(yr4_com.i), 
        w1  = rad2deg(apophis_com.w), 
        w2  = rad2deg(yr4_com.w),
        Om1 = rad2deg(apophis_com.Om), 
        Om2 = rad2deg(yr4_com.Om), 
        moid = NaN,
    ),
    (
        tag = "Earth-3I Atlas",
        rp1 = earth_com.q, 
        e1  = earth_com.e, 
        rp2 = atlas_3I_com.q, 
        e2  = atlas_3I_com.e, 
        i1  = rad2deg(earth_com.i), 
        i2  = rad2deg(atlas_3I_com.i), 
        w1  = rad2deg(earth_com.w), 
        w2  = rad2deg(atlas_3I_com.w),
        Om1 = rad2deg(earth_com.Om), 
        Om2 = rad2deg(atlas_3I_com.Om), 
        moid = NaN,
    ),
]

# orbit 1 for Wisniowsky tests
lutetia = (q = 2.036, e = 0.164, i = 0.0, w = 250.227, Om = 0.0)
for (i, row) in enumerate(eachrow(CSV.read(joinpath(@__DIR__, "orbit_tests.csv"), DataFrame)))
    push!(cases, (
            tag = "Wisniowsky $i",
            rp1 = lutetia.q,
            e1  = lutetia.e,
            rp2 = row["q_au"],
            e2  = row["e"],
            i1  = lutetia.i,
            i2  = row["i_deg"],
            w1  = lutetia.w,
            w2  = row["w_deg"],
            Om1 = lutetia.Om,
            Om2 = row["Om_deg"],
            moid = row["moid_au"]
        )
    )
end

""" debug build
"""
function moid_scan(
    ellipse1::Kepler.Ellipse{T, V}, 
    ellipse2::Kepler.Ellipse{T, V};
    sink = Any[],
    tol1::T  = 1e-4,
    tol2::T  = 1e-8, 
    max_iter = 100,
    ds::T    = 0.2,
    tol_cond = 1e3,
    ) where {T, V}
    # scanning pass, get initial guesses
    _, u_min, v_min = Kepler.moid_scan_meridional(ellipse1, ellipse2)

    # run levenberg-marquardt on each guess, take the lowest
    moid = Inf
    vf   = Inf
    uf   = Inf
    _f   = ((u, v),) -> sum(Kepler.dists_with_partials(u, ellipse1, v, ellipse2)[1] .^ 2)
    for (v, u) in zip(v_min, u_min)
        if isinf(v) || isinf(u)
            continue
        end

        push!(sink, [])
        push!(sink[end], (
                i = 0, 
                x = [u, v], 
                # dx = (du, dv), 
                # points = points,
                cost = _f([u, v]), 
                # mu = mu,
                # k = k,
            ))

        # use a simplex method to refine the initial guess
        points = MVector{3, SVector{2, T}}(
            SVector{2}(u + ds, v),
            SVector{2}(u, v + ds),
            SVector{2}(u - ds, v - ds),
        )
        vals = _f.(points)

        d0 = Inf
        x  = sum(points) ./ length(points)
        d  = _f(x)
        i  = 1
        push!(sink[end], (
            i = i, 
            x = [u, v], 
            # dx = (du, dv), 
            # points = points,
            cost = _f(x), 
            # mu = mu,
            # k = k,
        ))
        while abs(d0 - d) > tol1 && i < max_iter
            i += 1
            Kepler.simplex_step!(_f, points, vals)
            d0 = d
            x  = sum(points) ./ length(points)
            d  = _f(x)

            

            push!(sink[end], (
                i = i, 
                x = x, 
                # dx = (du, dv), 
                # points = points,
                cost = _f(x), 
                tag = "simplex",
                del_d = abs(d - d0),
                # mu = mu,
                # k = k,
            ))
        end
        
        # now run least squares until convergence (or worse)
        i    += 1
        u, v  = x
        dx, H = Kepler.dists_with_partials(u, ellipse1, v, ellipse2)
        Ht    = transpose(H)
        d     = dot(dx, dx)
        di    = Inf
        stepx = Inf
        while abs(d - d0) > tol2 && i < max_iter
            i += 1

            du, dv = inv(Ht*H)*Ht*dx
            dx, H  = Kepler.dists_with_partials(u + du, ellipse1, v + dv, ellipse2)
            Ht     = transpose(H)
            stepx  = sqrt(du^2 + dv^2)
            di     = dot(dx, dx)

            if di > d
                Kepler.simplex_step!(_f, points, vals)
                x  = sum(points) ./ length(points)
                di = _f(x)
                u, v = x

                push!(sink[end], (
                    i = i, 
                    x = (u, v), 
                    cost = _f(x), 
                    tag = "simplex",
                    del_d = abs(d - d0),
                ))

                dx, H = Kepler.dists_with_partials(u, ellipse1, v, ellipse2)
                Ht    = transpose(H)
            else
                v  = v + dv
                u  = u + du
                push!(sink[end], (
                    i = i, 
                    x = (u, v), 
                    cost = d, 
                    tag = "LS",
                    del_d = d0 - d,
                    # del_x = stepx,
                ))
            end

            # if cond(Ht*H) < tol_cond
            #     # reasonably well conditioned, continue normal equations solution
            #     du, dv = inv(Ht*H)*Ht*dx
            #     dx, H = Kepler.dists_with_partials(u + du, ellipse1, v + dv, ellipse2)
            #     Ht    = transpose(H)
            #     stepx = sqrt(du^2 + dv^2)
            #     di    = dot(dx, dx)



            #     v  = mod(v + dv, 2pi)
            #     u  = mod(u + du, 2pi)
            # else
            #     # take another simplex step
            #     Kepler.simplex_step!(_f, points, vals)
            #     x  = sum(points) ./ length(points)
            #     di = _f(x)

            #     u, v = x
            #     v  = mod(v, 2pi)
            #     u  = mod(u, 2pi)
            #     push!(sink[end], (
            #         i = i, 
            #         x = (u, v), 
            #         # dx = (du, dv), 
            #         # points = points,
            #         cost = _f(x), 
            #         tag = "simplex",
            #         del_d = abs(d - d0),
            #         # mu = mu,
            #         # k = k,
            #     ))
            #     # still need to compute matrix for next step
            #     dx, H = Kepler.dists_with_partials(u, ellipse1, v, ellipse2)
            #     Ht    = transpose(H)
            # end

            d0  = d
            d   = di
        end

        if d < moid
            moid = d
            vf   = v
            uf   = u
        end
    end

    return sqrt(moid), uf, vf
end


for (c, case) in enumerate(cases)
    orbit1 = Kepler.CometaryOrbit(case.rp1, case.e1, deg2rad(case.i1), deg2rad(case.Om1), deg2rad(case.w1), 0.0, 0.0, gm)
    orbit2 = Kepler.CometaryOrbit(case.rp2, case.e2, deg2rad(case.i2), deg2rad(case.Om2), deg2rad(case.w2), 0.0, 0.0, gm)

    ellipse1 = Kepler.Ellipse(orbit1)
    ellipse2 = Kepler.Ellipse(orbit2)

    # d, uf, vf = Kepler.moid_scan(ellipse1, ellipse2; mu = 100.0, tol = 1e-14)
    d, uf, vf = Kepler.moid_simplex(ellipse1, ellipse2; tol1 = 1e-5, tol2 = 1e-10)
    @printf "--- Case %d (%s) ---\n" c case.tag
    @printf "  MOID          = %12.3e au (%12.3e km)\n" d d*au
    @printf "  angle 1       = %11.6f deg (%11.6f deg) \n" rad2deg(uf) < 180 ? rad2deg(uf) : -(360 - rad2deg(uf)) rad2deg(uf)
    @printf "  angle 2       = %11.6f deg (%11.6f deg) \n" rad2deg(vf) < 180 ? rad2deg(vf) : -(360 - rad2deg(vf)) rad2deg(vf)
    @printf "  MOID - Ref    = %12.3e au\n" d - case.moid
    @printf "  Relative diff = %12.3e\n" (d - case.moid)/case.moid
    @printf "  rel err       = %12.3e\n" 2(d - case.moid)/(case.rp1 + case.rp2)
    println()
end

# 23 is VERY poorly behaved
# c = 5
for c in [5, 6, 7, 8]
    case = cases[c]
    orbit1 = Kepler.CometaryOrbit(case.rp1, case.e1, deg2rad(case.i1), deg2rad(case.Om1), deg2rad(case.w1), 0.0, 0.0, gm)
    orbit2 = Kepler.CometaryOrbit(case.rp2, case.e2, deg2rad(case.i2), deg2rad(case.Om2), deg2rad(case.w2), 0.0, 0.0, gm)

    ellipse1 = Kepler.Ellipse(orbit1)
    ellipse2 = Kepler.Ellipse(orbit2)

    u_grid = (0:0.001:1) .* 2pi
    v_grid = (0:0.001:1) .* 2pi
    M = zeros(Float64, length(u_grid), length(v_grid))
    @showprogress for ((i1, u), (i2, v)) in Iterators.product(enumerate(u_grid), enumerate(v_grid))
        dx, _ = Kepler.dists_with_partials(u, ellipse1, v, ellipse2)
        # M[i1, i2] = sqrt(dot(dx, dx))/2
        M[i1, i2] = norm(dx)
    end

    # f  = Figure()
    # ax = Axis(f[1, 1])
    # hm = contourf!(ax, rad2deg.(u_grid), rad2deg.(v_grid), M; levels = 50)
    # Colorbar(f[1, 2], hm)
    # f

    dmin, u_min, v_min = Kepler.moid_scan_meridional(ellipse1, ellipse2)
    # sink = []
    # d, uf, vf = Kepler.moid_scan(ellipse1, ellipse2; mu = 100.0, tol = 1e-14, nu = 0.0, k = 1.0)
    # case.moid - d

    sink = []
    d, uf, vf = moid_scan(ellipse1, ellipse2; sink = sink, tol1 = 1e-5, tol2 = 1e-12, max_iter = 200)
    d - case.moid

    # d, uf, vf = Kepler.moid_simplex(ellipse1, ellipse2; tol2 = 1e-10)
    # d - case.moid

    length.(sink)

    begin
    f  = Figure(size = (600, 500))
    ax = Axis(f[1, 1], xlabel = "f₁ (deg)", ylabel = "f₂ (deg)", title = case.tag,
        xtickalign = 1,
        ytickalign = 1,
    )

    cmap   = :viridis
    crange = (0, 5)

    contour!(ax, 
        rad2deg.(u_grid), 
        rad2deg.(v_grid), 
        M; 
        # levels = vcat(0:0.2:10, 10:30), 
        # levels = vcat(0:0.02:1, 1:0.1:4), 
        levels = vcat(0:0.05:0.2, 0.2:0.2:5),
        # levels = 50,
        # extendlow = :cyan, 
        # extendhigh = :magenta
        colormap = cmap,
        colorrange = crange,
    )
    Colorbar(f[1, 2]; label = "|distance| (au)", colormap = cmap, colorrange = crange)

    markers = [:circle, :utriangle, :rect, :star5]
    iter = 1
    sc = for (s, m, u, v) in collect(zip(sink, markers, u_min, v_min))
        if u == Inf || v == Inf
            continue
        end
        x = [mod(rad2deg(p.x[1]), 360) for p in s]
        y = [mod(rad2deg(p.x[2]), 360) for p in s]
        t = eachindex(x) ./ length(x)
        # scatterlines!(ax, x, y; markercolor = t, colormap = :plasma, marker = m)
        sc = scatter!(ax, x, y; color = t, colormap = :plasma, marker = m, markersize = 10.0, label = "Guess $iter")
        scatter!(ax, rad2deg(u), rad2deg(v); color = :red, markersize = 5.0)
        iter += 1
        sc
    end

    Colorbar(f[1, 3]; colormap = :plasma, label = "iteration (normalized)")

    scatter!(ax, rad2deg(uf), rad2deg(vf); color = :red, marker = :star5, label = "final MOID")
    xlims!(ax, 0, 360)
    ylims!(ax, 0, 360)
    ax.xticks = 0:30:360
    ax.yticks = 0:30:360
    axislegend(ax; position = :lt)
    # f
    save("/Users/samuelcornwall/school/courses/AE498-pd/HW/completed/HW1_plots/$(case.tag).png", f)
    end
end

# (uf, vf) = Kepler.simplex_solve(((u, v),) -> norm(Kepler.dists_with_partials(u, ellipse1, v, ellipse2)[1]), SVector{2}(u_min[1], v_min[1]); sink = sink)
# f = ((u, v),) -> norm(Kepler.dists_with_partials(u, ellipse1, v, ellipse2)[1])

# sink = []
# function nlls!(du, u, p)
#     du[:] .= Kepler.dists_with_partials(u[1], ellipse1, u[2], ellipse2)[1]
#     push!(sink, (u[1], u[2]))
# end

# u0   = [u_min[1], v_min[2]]
# prob = NonlinearSolve.NonlinearLeastSquaresProblem(
#     NonlinearSolve.NonlinearFunction(nlls!, resid_prototype = zeros(3)), u0)
# NonlinearSolve.solve(prob, NonlinearSolve.GaussNewton(), reltol = 1e-12, abstol = 1e-12)

# for x in sink[2:end]
#     scatter!(ax, x[1].value, x[2].value; color = :black)
# end


D = 2
x0 = [u_min[1], v_min[1]]
init_step = 0.1

f = ((u, v),) -> norm(Kepler.dists_with_partials(u, ellipse1, v, ellipse2)[1])

# --- START ---
begin
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "u (deg)", ylabel = "v (deg)")
hm = contourf!(ax, 
    rad2deg.(u_grid), 
    rad2deg.(v_grid), 
    M; 
    # levels = vcat(0:0.2:10, 10:30), 
    levels = 50, 
    # levels = 50,
    extendlow = :cyan, 
    extendhigh = :magenta
)
Colorbar(fig[1, 2], hm; label = "|distance| (au)")
fig
end

# --- First simplex
begin
points = [x0 for _ in 1:D+1]
for i in 1:D
    p = copy(x0)
    p[i] += init_step
    points[i] = p
    println(p)
end
points[end] = x0 .- init_step
vals = [f(p) for p in points]

_x = [rad2deg(x[1]) for x in points]
_y = [rad2deg(x[2]) for x in points]
scatter!(ax, _x, _y)
end

# --- Try a step ---
begin
simplex_step!(f, points, vals)
_x = [rad2deg(x[1]) for x in points]
_y = [rad2deg(x[2]) for x in points]
scatter!(ax, _x, _y)
end


d0 = f(x0)
d  = minimum(vals)
i  = 0

max_iter = 100
tol = 1e-8

while abs(d - d0) > tol && i < max_iter
# for i in 1:20
    i += 1
    println(i)

    println(points)
    _x = [rad2deg(x[1]) for x in points]
    _y = [rad2deg(x[2]) for x in points]
    lines!(ax, _x, _y)

    d0 = d
    d  = vals[j]
end


sink = []
d, uf, vf = Kepler.moid_scan(ellipse1, ellipse2; sink = sink, mu = 100.0, tol = 1e-14, nu = 0.0, k = 1.0)
case.moid - d

begin
f  = Figure()
ax = Axis(f[1, 1], xlabel = "u (deg)", ylabel = "v (deg)")
hm = contourf!(ax, 
    rad2deg.(u_grid), 
    rad2deg.(v_grid), 
    M; 
    # levels = vcat(0:0.2:10, 10:30), 
    levels = vcat(0:0.02:1, 1:0.1:4), 
    # levels = 50,
    extendlow = :cyan, 
    extendhigh = :magenta
)
Colorbar(f[1, 2], hm; label = "|distance| (au)")

markers = [:circle, :utriangle, :rect, :star]
for (s, m, u, v) in collect(zip(sink, markers, u_min, v_min))
    if u == Inf || v == Inf
        continue
    end
    x = [rad2deg(p.x[1]) for p in s]
    y = [rad2deg(p.x[2]) for p in s]
    t = eachindex(x) ./ length(x)
    # scatterlines!(ax, x, y; markercolor = t, colormap = :plasma, marker = m)
    scatterlines!(ax, x, y; markercolor = t, colormap = :plasma, marker = m)
    scatter!(ax, rad2deg(u), rad2deg(v); color = :red, markersize = 4.0)
end

scatter!(ax, rad2deg(uf), rad2deg(vf); color = :red, marker = :star5)
xlims!(ax, 0, 360)
ylims!(ax, 0, 360)
ax.xticks = 0:30:360
ax.yticks = 0:30:360
f
end