# %%
using Pkg; Pkg.activate()
using Kepler
using StaticArrays
using LinearAlgebra
using Printf
using SPICE
using Logging
using Roots
using BenchmarkTools
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

# pos = [-1.8334693215194804, -0.7424496266825248, -0.2952266528773983]
# vel = [0.004858030562369436, -0.010253429035065403, -0.004384412542682951]
# dt = 3.0695209246163175
# gm = 0.0002959122082326087

pos = [-2.605657, -1.4783354, 0.15835176]
vel = [0.004609774, -0.008403798, -0.0026038089]
dt  = -3.163144998718053
gm  = 0.0002959122082326087

  #   pos = Float32[-1.3784086, -2.4442484, -1.0609612]
#   vel = Float32[0.008740334, -0.0046784547, -0.0005773584]
#   dt  = 1.171095000114292
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.9948885, -1.8434341, -1.2736433]
#   vel = Float32[0.0073577547, -0.006099169, -0.0026961644]
#   dt  = -5.954300998710096
#   gm  = 0.0002959122082326087
#   pos = Float32[0.70663416, -2.758167, -0.9450836]
#   vel = Float32[0.009562632, 0.001755457, 0.0020267316]
#   dt  = -3.8347659991122782
#   gm  = 0.0002959122082326087
#   pos = Float32[-2.548434, -1.3640729, -0.80299234]
#   vel = Float32[0.005233072, -0.006991323, -0.0047299895]
#   dt  = -1.9699699995107949
#   gm  = 0.0002959122082326087
#   pos = Float32[-2.7282562, -1.2332435, -0.18902071]
#   vel = Float32[0.0040133526, -0.008320237, -0.003642466]
#   dt  = -2.059900999534875
#   gm  = 0.0002959122082326087
#   pos = Float32[-0.89066714, -2.6600819, -1.0633336]
#   vel = Float32[0.009474554, -0.0025725586, -0.0015003541]
#   dt  = 6.139203998725861
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.5921421, -2.3637965, -0.9367762]
#   vel = Float32[0.008296401, -0.0042097615, -0.0034778588]
#   dt  = 1.179552000015974
#   gm  = 0.0002959122082326087
#   pos = Float32[-2.4770048, -1.623562, -0.47800973]
#   vel = Float32[0.0055334233, -0.008206098, -0.00080279325]
#   dt  = -1.9189719995483756
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.3575367, -2.5892, -0.6731543]
#   vel = Float32[0.008848818, -0.004448635, -0.00073412963]
#   dt  = 3.1200919994153082
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.4613756, -2.561916, -0.54860514]
#   vel = Float32[0.00862502, -0.0049213725, 6.777164f-6]
#   dt  = 1.1746660000644624
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.8186375, -1.9565544, -1.3654505]
#   vel = Float32[0.007881786, -0.005304774, -0.0028969562]
#   dt  = -5.952940998598933
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.7796177, -2.1878915, -1.0227859]
#   vel = Float32[0.007855743, -0.0044626654, -0.0041225282]
#   dt  = 5.1331169991753995
#   gm  = 0.0002959122082326087
#   pos = Float32[-2.8041542, -1.0645479, 0.058787927]
#   vel = Float32[0.0032609692, -0.008773213, -0.0033213517]
#   dt  = -3.157348998822272
#   gm  = 0.0002959122082326087
#   pos = Float32[-2.5528412, -1.2610574, -0.944847]
#   vel = Float32[0.005042023, -0.00806005, -0.0028660588]
#   dt  = -2.0183119997382164
#   gm  = 0.0002959122082326087
#   pos = Float32[-1.0205196, -2.5464525, -1.2141327]
#   vel = Float32[0.0089850165, -0.0017688923, -0.003842312]
#   dt  = 6.12485099863261
#   gm  = 0.0002959122082326087

pos = [-0.0005851491971243794, -2.964878053101466, -1.155167716817128]
vel = [-9.165459464372124e-5, 0.04049187111210227, 0.015647875649757625]
dt = 9.96469299821183
gm = 0.0002959122082326087

# test and compare to spice
posf, velf = Kepler.propagate(pos, vel, dt, gm)
statef = SPICE.prop2b(gm, [pos..., vel...], dt)
posf2 = statef[1:3]
velf2 = statef[4:6]

posf - posf2
velf - velf2

# becnhmark
p0 = SVector{3}(pos)
v0 = SVector{3}(vel)

@benchmark Kepler.propagate($p0, $v0, $dt, $gm)

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