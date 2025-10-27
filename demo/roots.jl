using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Roots
using Logging
using CairoMakie

# test orbit
pos = [1.0, 0.0, 0.0]
# vel = [0.0, sqrt(2), 0.0] # parabolic
vel = [0.0, 1.5, 0.0] # hyperbolic
# vel = 
# gm  = 1.0 # assume normalized
dt = pi/2

# get the initial guess (bracket)
dr0 = dot(vel, pos)
a   = 2.0 - dot(vel, vel) # 1/semimajor axis! alpha = beta in these coordinates

Hvec = cross(pos, vel)
Evec = cross(vel, Hvec) - pos
e    = norm(Evec)

x = dt*a/(1 - e)
x_br = 0.0
y_br = -dt

# _, c1, c2, c3 = Kepler.stumpff(a*x^2)
# y = x*c1 + dr0*c2*x^2 + c3*x^3 - dt
g0, g1, g2, g3 = Kepler.stumpff(a, x)
y = g1 + dr0*g2 + g3 - dt

while sign(y) == sign(y_br)
    # take forward steps, maintaining the size of our interval
    # we could use some newton steps here, and would probably work better
    # but I'm lazy, and this is exceedingly likely to require at most 1 step
    x_br = x
    y_br = y
    x    *= 2
    g0, g1, g2, g3 = Kepler.stumpff(a, x)
    y    = g1 + dr0*g2 + g3 - dt
end

# we now have a bracket, and can find the root.
# chain a couple of closures together to get the right args to the right places
# TODO: implement a bracketed root finder using higher derivatives
bracket = (x_br, x)

# run a known root finder
_f1 = _x -> (_x, Kepler.stumpff(a*_x^2)...)
_f2 = ((_x, _, _c1, _c2, _c3),) -> _x*_c1 + dr0*_c2*_x^2 + _c3*_x^3
x_actual = find_zero(_x -> _f2(_f1(_x)) - dt, bracket, A42())

# try the new rootfinder
function univ_kepler_2(x, a, dt)
    # c0, c1, c2, c3 = Kepler.stumpff(a*x^2)
    g0, g1, g2, g3 = Kepler.stumpff(a, x) 
    f  = g1 + dr0*g2 + g3 - dt
    f1 = g0 + dr0*g1 + g2
    return f, f1
end

x_l = x_br
x_h = x
f_l, f1_l = univ_kepler_2(x_l, a, dt)
f_h, f1_h = univ_kepler_2(x_h, a, dt)

x_steps = Float64[]
f_steps = Float64[]
dx = Inf
fi = min(f_l, f_h)
# i = 0

# bootstrap with a 2 point step
# we can use higher derivatives if available
i = 1
x3 = Kepler.lmm12_step(x_l, x_h, f_l, f_h, f1_l, f1_h)
f3, f13, = univ_kepler_2(xi, a, dt)

# continue with 3 points
while abs(fi) > 0 && dx > 1e-16 && i < 100
    i += 1
    # xi = Kepler.lmm12_step(x_l, x_h, f_l, f_h, f1_l, f1_h)
    xi = Kepler.lmm12_step(x_l, x_3, x_h, f_l, f_3, f_h, f1_l, f13, f1_h)
    fi, f1i, = univ_kepler_2(xi, a, dt)

    push!(x_steps, xi)
    push!(f_steps, fi)

    # find the new bracket
    # from the LMM paper, past 3 points is probably not worth it.
    # so we retain the 3 most tightly bracketing the root
    if fi == 0
        # dx = 
        continue
    elseif sign(fi) == sign(f_l)
        dx   = abs(x_l - xi)
        x_l  = xi
        f_l  = fi
        f1_l = f1i
    elseif sign(fi) == sign(f_h)
        dx   = abs(x_h - xi)
        x_h  = xi
        f_h  = fi
        f1_h = f1i
    else
        throw("AHHHHHH")
    end
    # println()
    @printf "x: %.6e f: %.6e dx: %.6e\n" xi fi dx
end

begin
    fig  = Figure(size = (500, 800))
    ax   = Axis(fig[1, 1], xlabel = "x", ylabel = "dt - dt0")
    ax2  = Axis(fig[2, 1], xlabel = "x", ylabel = "(dt - dt0)/ds")

    xi  = -0.5:0.01:2
    yi  = [univ_kepler_2(xii, a, dt)[1] for xii in xi]
    f1  = [univ_kepler_2(xii, a, dt)[2] for xii in xi]
    # y1i = [univ_kepler_2(xii, a, dt)[2] for xii in xi]
    # zi  = -2.0:0.01:0.5
    # y2i = b0 .+ (b1 .* zi) .+ (b2 .* zi .^ 2) .+ (b3 .* zi .^ 3)

    lines!(ax, xi, yi)
    lines!(ax2, xi, f1)
    # lines!(ax, y2i, zi)
    scatter!(ax, x_br, y_br)
    scatter!(ax, x, y)
    # scatter!(ax, x2, f_2; marker = :utriangle, markersize = 20)
    # scatter!(ax, x2, f_2; marker = :utriangle, markersize = 20)
    for (xi, fi) in zip(x_steps, f_steps)
        scatter!(ax, xi, fi; marker = :utriangle, markersize = 20)
    end
    # scatter!(ax, x_actual, 0.0)

    # xlims!(ax, 1.4, 1.45)

    fig
end
