using Kepler
using LinearAlgebra

pos = [0.9921454991096386, 0.17725124417460564, 0.0]
vel = [-0.0021389710259218375, 0.024134905635018196, 0.0]
dt = 3.653014713476504
gm = 0.00029584

r0 = norm(pos)
ts = sqrt(r0^3 / abs(gm))
T = dt / ts # non-dim time of flight
b = 2.0 - dot(vel, vel)*r0/gm   # non-dim negative energy
k = dot(vel, pos) / sqrt(r0*gm) # non-dim r . v

P  = 2pi/sqrt(abs(b)^3)
dT = T + P/2 - 2k/b
n  = 0
if b > 0
    n = floor(Int, dT/P)
    T -= n*P
end

T
T/P

u0 = 0.0
t0, df0 = Kepler.tof_nondim(u0, b, k)
f0 = t0 - T

s1 = T
# # s1 = sign(T)*abs(T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4))
# s1 = if abs(T) < abs(P)/4 || b > 1e-6
#     # sign(T)*abs(T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4))
#     T - k/2*(T^2) - (1 - b - 3k^2)/6*(T^3) + k*(10 - 9*b - 15*k)/24*(T^4)
# else
#     if b < -1e-6
#         1.01sign(T)*log(
#             (-2*b*T)/(k + sign(T)*(1 - b)/sqrt(-b))
#         )/sqrt(-b)
#         # -sign(T)*log(abs(2b*T/(k - (1 - b)/sqrt(-b))))/sqrt(-b)

#         # sign(T)*log(
#         #     abs(b*T*(
#         #         1/(k + sign(T)*(1 - b)/sqrt(-b)) +
#         #         1/(k - sign(T)*(1 - b)/sqrt(-b)) +
#         #         1
#         #     ))
#         # )/sqrt(-b)

#     else
#         hvec = cross(pos, vel)
#         p = dot(hvec, hvec)
#         z = acot(3*T*sqrt(1/p^3))/2
#         w = atan(cbrt(tan(z)))
#         2*sqrt(p)*cot(2w)
#     end
# end


u1 = if b > 0
    tan(sqrt(b)*abs(s1)/4)/sqrt(b)
elseif b < 0
    tanh(sqrt(-b)*abs(s1)/4)/sqrt(-b)
elseif b == 0
    abs(s1)/4
else
    throw("invalid b")
end
u1 = sign(T)*min(abs(u1), 1/sqrt(abs(b))) # we can bound u just in case the guess overshoots wildly


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


ui = if u0 == u1
    u0
elseif f1 == 0
    u1
else
    Kepler.flmsm1_step((x = u0, y = f0, dy = df0), (x = u1, y = f1, dy = df1))
end

if !(min(u0, u1) <= ui <= max(u0, u1))
    ui = (u0 + u1)/2
end


i = 0
du2 = Inf
du1 = abs(u0 - u1)
while ui != u0 && ui != u1
    i += 1
    if i == 2000
        throw("max iterations")
    end
    ti, dfi = Kepler.tof_nondim(ui, b, k)
    fi = ti - T
    # shift appropriate end point
    if fi == 0
        break
    elseif sign(fi) == sign(f0)
        u0, f0, df0 = (ui, fi, dfi)
    elseif sign(fi) == sign(f1)
        u1, f1, df1 = (ui, fi, dfi)
    else
        throw("invalid sign detected $(fi)")
    end
    # how much is the interval shrinking
    du0 = abs(u0 - u1)
    
    # compute next point
    # ui = Kepler.flmsm1_step((x = u0, y = f0, dy = df0), (x = u1, y = f1, dy = df1))
    dfi = (f1 - f0)/(u1 - u0)
    ui -= fi/dfi

    if !(min(u0, u1) <= ui <= max(u0, u1))
        ui = (u0 + u1)/2
        # ui = u0 + (u1 - u0)*rand()/2
        println("bisect (missed)")
    elseif du0 > du2/2
        ui = (u0 + u1)/2
        # ui = u0 + (u1 - u0)*rand()/2
        println("bisect (slow)")
    else
        println("interpolate")
    end

    du2 = du1
    du1 = du0
end
u0, u1
f0, f1

ui
ti, dfi = Kepler.tof_nondim(ui, b, k)
fi = ti - T

dfi = (f1 - f0)/(u1 - u0)
fi / dfi


i


# secant u    = 0.015455209558264502
# interpolate = 0.015455209558264502