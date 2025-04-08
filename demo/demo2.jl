# %%
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
# %%

# a failure case encountered using heliolinc. hyperoblic, large |s|, sinh(sqrt(|z|)) overflows
# turns out it was a negative time issue (and the overflowed initial guess). chooseing a better guess and reversing time fixed it
# %%
# pos = [2.5, 0.0, 0.0]
# vel = [-0.0185, 4.8805672928708536e-5, 0.0]
pos = [2.5, 0.0, 0.0]
vel = [-0.015, 0.003378725763145275, 0.0]
dt = 1.234750019852072
# 
# 0.0002959122082326087
# dt  = -11.212502314709127
gm  = 0.01720209894846^2
# %%

posf, velf = Kepler.solve(pos, vel, dt, gm)
vcat(posf, velf) - SPICE.prop2b(gm, [pos..., vel...], dt)



r  = norm(pos)
dr = dot(vel, pos) / r
v2 = dot(vel, vel)
alpha = 2gm/r - v2
hvec = cross(pos, vel)
evec = cross(vel, hvec)/gm - pos/r
e = norm(evec)

# s0 = dt*alpha/(gm*(1-e))

# trying new hyperbolic guess 
k  = 1.8
dM = sqrt(-gm*alpha^3)*dt
CH = 1 - r*alpha
SH = dot(pos, vel)*sqrt(-alpha/gm)
dF = dM > 0 ? log((2dM + k*e)/(CH + SH)) : log((-2dM + k*e)/(CH - SH)) 
s0 = dF / sqrt(-alpha)

Kepler.universal_kepler(0.0, alpha, r, dr, gm, 0.0)
Kepler.universal_kepler(s0, alpha, r, dr, gm, 0.0)

z = alpha*s0^2 # very big

n = 0
while abs(z) > 10.0
    z /= 4
    n += 1
end

z, n

sinh2 = sinh(0.5sqrt(-z))
cosh2 = cosh(0.5sqrt(-z))
c1 = 2sinh2*cosh2/sqrt(-z)
c2 = -2sinh2^2/z
c3 = (1 - c1)/z
c0 = 1 - z*c2
c0, c1, c2, c3

while n > 0
    c32 = (c2 + c0*c3)/4
    c22 = (c1^2)/2
    c12 = c0*c1
    c02 = 2(c0^2) - 1
    n -= 1
    c0 = c02
    c1 = c12
    c2 = c22
    c3 = c32
    println("$(n) $(c0) $(c1) $(c2) $(c3)")
end

