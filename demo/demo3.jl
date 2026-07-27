# Similarly to the conservation plots, we scan the orbit space and compare to the spice results
using SPICE
using Kepler
using LinearAlgebra
using Printf
using StaticArrays
using GLMakie

gm  = 0.0172^2
pos = [1.0, 0.0, 0.0]
vel = [0.0, sqrt(gm)*1.1, 0.0]

a = 1/(2 - dot(vel, vel)/gm) # r = 1

hvec = cross(pos, vel)
evec = cross(vel, hvec)/gm - normalize(pos)

e = norm(evec)
@info @sprintf "a = %.6e e = %.6f" a e

n   = sqrt(gm)
P   = 2pi*sqrt(abs(a^3/gm))

x = Float64[]
y = Float64[]
z = Float64[]
t = (0:0.01:1) .* P

for dt in t
    posf, _ = Kepler.propagate_sheppard(pos, vel, dt, gm)
    # push!(traj, SVector{3}(posf))
    push!(x, posf[1])
    push!(y, posf[2])
    push!(z, posf[3])
end

f  = Figure()
ax = Axis3(f[1, 1])
scatter!(ax, x, y, z; color = t)

l = 2.0
for lims! in (xlims!, ylims!, zlims!)
    lims!(ax, -l, l)
end
f

dt = 0.75*P
# debug
r0 = norm(pos)
ts = sqrt(r0^3 / abs(gm))
T = dt / ts # non-dim time of flight
b = 2.0 - dot(vel, vel)*r0/gm   # non-dim negative energy
k = dot(vel, pos) / sqrt(r0*gm) # non-dim r . v

P  = 2pi/sqrt(abs(b)^3)
dT = T + P/2 - 2k/b
n  = if b > 0
    floor(Int, dT/P)
else
    0
end

T -= n*P

T/P

u0 = 0.0
t0, df0 = Kepler.tof_nondim(u0, b, k)
f0 = t0 - T

s1 = T
u1 = if b >= 0
    tan(sqrt(b)*abs(s1)/4)/sqrt(b)
else
    tanh(sqrt(-b)*abs(s1)/4)/sqrt(-b)
end

u1, 1/sqrt(abs(b))

u1 = sign(T)*min(abs(u1), 1/sqrt(abs(b)))

t1, df1 = Kepler.tof_nondim(u1, b, k)
f1 = t1 - T

f0, f1

i = 0
while i < 10 && sign(f0)*sign(f1) > 0 #|| !isfinite(f0) || !isfinite(f1))
    i += 1
    if isfinite(f1)
        # initial guess fell short, we can shift
        (t0, f0, df0) = (t1, f1, df1)
        u1 = (u1 - u0) + u1
        # if abs(u1) > 1/sqrt(abs(b))
        #     u1 = sign(T)/sqrt(abs(b))
        # end
    elseif isfinite(f0) 
        # overshot so far that the guess overflows, have to scale back
        u1 = (u0 + u1)/2
    else
        # both are infinite, we cannot find the bracket
        throw(DomainError(dt, "The time span exceeds the provided precision for this orbit. Consider normalizing to a different scale or taking intermediate steps."))
    end
    t1, df1 = Kepler.tof_nondim(u1, b, k)
    f1 = t1 - T
    println((u0, u1, f0, f1))
end

u0, u1
f0, f1

exit()

vr_range = 10 .^ (-5:0.2:-4)
# vr = vcat(-reverse(vr_range), 0.0, vr_range)
vr = 0.0:0.0
vt = 10 .^ (-4:0.1:-1)
# vt = 0.0:0.0

dt = 1e-3*T

for (vri, vti) in Iterators.product(vr, vt)
    vel = [vri, vti, 0.0]

    a = 1/(2 - dot(vel, vel)/gm) # r = 1

    hvec = cross(pos, vel)
    evec = cross(vel, hvec)/gm - normalize(pos)

    e = norm(evec)
    @info @sprintf "a = %.6e e = %.6f" a e

    try
        # pfk, vfk = Kepler.propagate(pos, vel, dt, gm)
        pfk, vfk = Kepler.propagate_sheppard(pos, vel, dt, gm)
        pvfs = SPICE.prop2b(gm, [pos..., vel...], dt)
        @info @sprintf "  position error: %.6e %.6e %.6e (%.6e)" (pvfs[1:3] .- pfk)... norm(pvfs[1:3] .- pfk)
    catch _
        throw("vr = $vri, vt = $vti")
    end

end


for (rp, e) in Iterators.product(
    10 .^ (-3:0.1:2),
    0:0.01:1
)
    a = rp/(1 - e)
    @info @sprintf "a = %.3e e = %.3f" a e
    n   = sqrt(gm)
    T   = 2pi/n

    v  = sqrt(gm*(2 - 1/a))
    vt = sqrt(gm*abs(a)*(1 - e^2))
    vr = sqrt(v^2 - vt^2)

    vel = [vr, vt, 0.0]
        #    n += 1
    for dt in 10 .^ (-3:0.1:0)
        pfk, vfk = Kepler.propagate(pos, vel, dt*T, gm)
        pvfs = SPICE.prop2b(gm, [pos..., vel...], dt*T)
        @info @sprintf "position error: %.6e %.6e %.6e" (pvfs[1:3] .- pfk)...
    end
end