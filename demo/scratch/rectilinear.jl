using Kepler
using CairoMakie

using Roots

b  = 2.0
s0 = 0.0

function tfunc(x, b, s0)
    _, U1, U2, U3 = Kepler.universal03(b, x)
    return U1 + s0*U2 + U3
end

function rfunc(x, b, s0)
    U0, U1, U2, _ = Kepler.universal03(b, x)
    return U0 + s0*U1 + U2
end


# x0 = Roots.findzero()
x1 = 0.0
x2 = 4.0
x_min = 2.0
r_min = rfunc(x_min, b, s0)

r1 = 1.0
r2 = rfunc(x2, b, s0)

dr = 
# while r_min > 1e-14
#     if 
# end



x = 0:0.001:5
y = rfunc.(x, b, s0)

f  = Figure()
ax = Axis(f[1, 1], xlabel = "s", ylabel = "r")
lines!(ax, x, y)
f