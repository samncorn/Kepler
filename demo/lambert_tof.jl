using Kepler
using StaticArrays
using LinearAlgebra
using CairoMakie


# phi = -4pi
r1 = 3.0
# nondimensionlaize to r1 = 1.0 or r2 = 1.0
r2 = 1.0
dv = pi/26
cosdv = cos(dv)
# intermediate variable
A = sqrt(r1*r2*(1 + cosdv))

n = 2
z = -100:0.1:400
tof_long  = [NaN for _ in z]
tof_short = [NaN for _ in z]

for (i, zi) in enumerate(z)
    try 
        t = Kepler.universal_lambert(zi, A, r1, r2, 1.0)
        tof_short[i] = t
    catch e
        
    end

    try 
        t = Kepler.universal_lambert(zi, -A, r1, r2, 1.0)
        tof_long[i] = t
    catch e
        # continue
    end
end

f  = Figure()
ax = Axis(f[1, 1])

lines!(ax, z, tof_short; label = "short")
lines!(ax, z, tof_long; label = "long")

ylims!(ax, -1, 400)
axislegend(ax)

ax.xticks = -100:50:400
ax.yticks = -100:20:400

f
