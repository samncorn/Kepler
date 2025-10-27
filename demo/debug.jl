# %%
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
using Roots
# %%

debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger)

# test from vallado
# T = Float64
# pos = T.([1131.340, -2282.343, 6672.423])
# vel = T.([-5.64305, 4.30333, 2.42879])
# dt  = T(40.0*60.0)
# gm  = T(398600.4415)# per vallado

# case that fails in the wild [FIXED]
# pos = [1.25, 0.0, 0.0]
# vel = [-0.021759125285088412, -0.015386025041773686, -0.0]
# dt  = 4.0
# gm  = 0.0002959122082326087

# parabolic orbit that failed with elliptic upper bound [FIXED]
# pos = [1.25, 0.0, 0.0]
# vel = [-0.010879562642544206, -0.018843955261014882, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# # near rectilinear, consider as invalid input
# pos = [2.5, 0.0, 0.0]
# vel = [0.007693012520886843, 1e-11, -0.0]
# dt  = 8.0
# gm  = 0.0002959122082326087

# failed in the wild
# runs fine, probably forgot to update to latest version 
# pos = [0.4186942098016838, 0.24637616481752644, 1.3672151662793137]
# vel = [0.023317918469691952, -0.014641183526584653, 0.0008003774012850409]
# dt  = -5.998977942803302
# gm  = 0.0002959122082326087

pos = [-1.8334693215194804, -0.7424496266825248, -0.2952266528773983]
vel = [0.004858030562369436, -0.010253429035065403, -0.004384412542682951]
dt = 3.0695209246163175
gm = 0.0002959122082326087

# test and compare to spice
posf, velf = Kepler.propagate(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)
posf2 = statef[1:3]
velf2 = statef[4:6]

posf - posf2
velf - velf2

# # break down internals
function stumpff(z)
    if z > 1e-6
        c0 = cos(sqrt(z))
        c1 = sin(sqrt(z))/sqrt(z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        return c0, c1, c2, c3
    elseif z < -1e-6
        c0 = cosh(sqrt(-z))
        c1 = sinh(sqrt(-z))/sqrt(-z)
        c2 = (1 - c0)/z
        c3 = (1 - c1)/z
        return c0, c1, c2, c3
    else
        # z very small. evaluate the series to 6 terms 
        # error is O(x^7) < 1e-21, well within floating point tolerances
        c2 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/182)/132)/90)/56)/30)/12)/2
        c3 = (1-z*(1-z*(1-z*(1-z*(1-z*(1-z/210)/156)/110)/72)/42)/20)/6
        c0 = 1 - z*c2
        c1 = 1 - z*c3
        return c0, c1, c2, c3
    end
end

function universal_kepler(y, l, k1, k2, k3)
    _, c1, c2, c3 = stumpff(l*y^2)
    L = y*(k1*c1 + y*(k2*c2 + y*k3*c3))
    return L
end

hvec = cross(pos, vel)
h2   = dot(hvec, hvec)
r0   = norm(pos)
dr0  = dot(vel, pos)/r0
v20  = dot(vel, vel)
a    = 2gm/r0 - v20
e    = sqrt(1 - a*h2/gm^2)
q    = h2/(gm*(1 + e))
l    = (1-e)/(1+e)

nu = sqrt(gm*(1 + e)/(q^3))
L  = nu*dt
k1 = r0/q
k2 = r0*dr0/(nu*(q^2))
k3 = gm/((nu^2)*(q^3))

bracket = (0.0, L)
x0 = sqrt(q/(1+e))*L
z0 = l*L^2 

universal_kepler(1.1L, l, k1, k2, k3) - L

y = find_zero(_y -> universal_kepler(_y, l, k1, k2, k3) - L, L, A42())
universal_kepler(y, l, k1, k2, k3) - L

# 
x = sqrt(q/(1+e))*y
l*y^2