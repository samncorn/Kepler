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
using ForwardDiff

import Test: @inferred
using Profile
# %%

# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debug_logger)

# test from vallado
T = Float64
pos = T.([1131.340, -2282.343, 6672.423])
vel = T.([-5.64305, 4.30333, 2.42879])
dt  = T(40.0*60.0)
gm  = T(398600.4415) # per vallado

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

# pos = [-2.605657, -1.4783354, 0.15835176]
# vel = [0.004609774, -0.008403798, -0.0026038089]
# dt  = -3.163144998718053
# gm  = 0.0002959122082326087

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

# [FIXED] near rectilinear, need to devise a solution 
# pos = [-0.0005851491971243794, -2.964878053101466, -1.155167716817128]
# vel = [-9.165459464372124e-5, 0.04049187111210227, 0.015647875649757625]
# dt = 9.96469299821183
# gm = 0.0002959122082326087

# [FIXED?] tring an actiual rectilinear case
# SPICE cannot handle rectilinear cases, but the universal variables equations are still valid
# pos = [-0.0005851491971243794, -2.964878053101466, -1.155167716817128]
# vel = [0.0, 0.0, 0.0]
# dt = 9.96469299821183
# gm = 0.0002959122082326087

# found in the wild
# very hyperbolic, initial guess overflows
# pos = [2.0018040983461143, 4.310971608104017, -1.4324077743346348]
# vel = [340.5788770104138, 455.2580303275891, -243.9375394403128]
# dt = 9.582662548925576
# gm = 0.0002959122082326087

# pos = [2.2275391726284246, 57.18786513542615, 655.9177224430462]
# vel = [0.999466534019431, 30.149863863535487, -878.3833478838895]
# dt  = 630.2967702062505

# pos = [109.62369242851923, 67.80081611193758, 49.75433855353975]
# vel = [-139.237922907874, -87.24739293809648, -62.59488732836259]
# dt = 5.949454855270001
# gm = 0.0002959122082326087

# pos = [40401.932642594336, -47224.957914188126, -33526.83950111518]
# vel = [42486.22331901802, -4893.683561918934, 22635.55566101651]
# dt  = 6.121242998633534
# # dt = -526779.8806769907
# gm  = 0.0002959122082326087

# pos = [41030.87830706519, -2349.8572376691723, 55434.04184220392]
# vel = [-7.565476082816816e6, 541902.5735269253, 2.8492970861261785e6]
# dt = 5.95204599853605
# gm = 0.0002959122082326087

# pos = [261.25387200460534, -211.11168024875982, 71.14807690957349]
# vel = [-11138.034281711645, -21791.64218327744, -9196.569919801788]
# dt = 71417.14474731834
# gm = 0.0002959122082326087

# pos = [-1.0869055514657402, -1.5771341941111423, -0.5908430891461786]
# vel = [0.024192638271511267, 0.014795405462042815, 0.004641678916793392]
# dt  = 1.178212999831885
# gm  = 0.0002959122082326087

# pos = [3.48963168414826, 6.823118245748422, 3.533810926367658]
# vel = [-0.007296823907528163, -0.12053259794490012, -0.051404238147262506]
# dt = 1.1817809999920428
# gm = 0.0002959122082326087

# benchmark
function test()
    _pos = SVector{3}(1131.340, -2282.343, 6672.423)
    _vel = SVector{3}(-5.64305, 4.30333, 2.42879)
    _dt  = 40.0*60.0
    _gm  = 398600.4415

    pos, vel = Kepler.propagate(_pos, _vel, _dt, _gm)
    # return pos, vel
end

@code_warntype test()
@btime test()
@allocated test()
# @profile test()
# Profile.print()

# posf, velf, dxdx, dxdv, dvdx, dvdv = Kepler.propagate_with_partials(pos, vel, dt, gm)
posf, velf = Kepler.propagate(pos, vel, dt, gm)

# check state against spice
statef = SPICE.prop2b(gm, [pos..., vel...], dt)
posf2 = statef[1:3]
velf2 = statef[4:6]

posf - posf2
velf - velf2

# check partials
dxdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[1], pos)
dxdx .- dxdx_auto

dxdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[1], vel)
dxdv .- dxdv_auto

dvdx_auto = ForwardDiff.jacobian(x -> Kepler.propagate(x, vel, dt, gm)[2], pos)
dvdx .- dvdx_auto

dvdv_auto = ForwardDiff.jacobian(x -> Kepler.propagate(pos, x, dt, gm)[2], vel)
dvdv .- dvdv_auto

# check orbital elements
q, e, i, Om, w, tp = Kepler.cometary(pos, vel, dt, gm)

if dt < 0
    dt  = -dt
    vel = -vel
end

# unwind the function to find problem
r0 = norm(pos)
s0 = dot(pos, vel)
b  = 2gm/r0 - dot(vel, vel)

# try to bracket the root
# dt/dx = r >= 0, and monotonic, so we can step until we find a bracket
# bracket = (0.0, dt/r)
xl = 0.0
yl = -dt
rl = r0

xh = if abs(gm*b) < 1e-5
    # parabolic (Vallado)
    h = cross(pos, vel)
    p = dot(h, h)/gm
    s = acot(3*sqrt(gm/p^3)*dt)/2
    w = atan(cbrt(tan(s)))
    sqrt(p)*2*cot(2w)/sqrt(gm)
elseif gm*b < 0
    # hyperbolic (Vallado)
    a = gm/b
    sqrt(-a)*log(-2gm*dt/(a*(s0+sqrt(-gm*a)*(1 - r0/a))))/sqrt(gm)
elseif gm*b > 0
    # elliptic
    dt/r0
end

dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)
yh = dth - dt

# xh = (xl + xh)/2
# Kepler.Kepler.stumpff(b*xh^2)
# dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)

xl, xh
yl, yh

b*xh^2

i = 0
while i < 100 && (sign(yl) == sign(yh) || isinf(dth) || isnan(dth))
    i += 1
    println(i)
    if isinf(dth) || isnan(dth) #|| isnan(rh)
        println("bisecting")
        xh = (xl + xh)/2
        dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)
        yh = dth - dt
        if xl == xh 
            throw("dt exceeds the computable range of values")
        end
    elseif sign(yl) == sign(yh) && !isnan(rh)
        println("shifting bracket")
        xl  = xh
        # xh += (dt - dth)/rh
        xh *= 2
        yl  = yh
        rl  = rh
        dth, rh = Kepler.universal_kepler2(xh, b, r0, s0, gm)
        yh  = dth - dt
    end
end

x = 0.5*(xh + xl)
i = 0
tol = 0.0
while abs(xh - xl) > tol && i < 10000
    # x = 0.5(xl + xh)
    i += 1
    x = Kepler.lmm12_step(xl, xh, yl, yh, rl, rh)
    if x == xl || x == xh
        break
    end
    y, r = Kepler.universal_kepler2(x, b, r0, s0, gm)
    y -= dt
    if sign(y) == sign(yl)
        xl = x
        yl = y
        rl = r
    elseif sign(y) == sign(yh)
        xh = x
        yh = y
        rh = r
    elseif sign(y) == 0
        break
    end
end
x = 0.5(xl + xh)
i


# (xl, xh)
# xh - xl
# Kepler.universal_kepler(xl, b, r0, s0, gm) - dt
# @benchmark Kepler.universal_kepler($xh, $b, $r0, $s0, $gm) - $dt

# @benchmark x = find_zero(_x -> Kepler.universal_kepler(_x, $b, $r0, $s0, $gm) - $dt, ($xl, $xh), A42())
# @benchmark x = (try
#     method = $A42()
#     # method = Bisection()
#     find_zero(_x -> Kepler.universal_kepler(_x, $b, $r0, $s0, $gm) - $dt, ($xl, $xh), method)
# catch _
#     # throw((pos = $pos, $vel = vel, dt = dt, gm = gm, bracket = (xl, xh)))
#     throw("HELL")
# end)

z = b*x^2
_, c1, c2, c3, c4, c5 = Kepler.stumpff5(b*x^2)
f    = 1 - (gm/r0)*(x^2)*c2
g    = dt - gm*(x^3)*c3
posf = f*pos + g*vel
rf   = norm(posf)
df   = -(gm/(rf*r0))*x*c1
dg   = 1 - (gm/rf)*(x^2)*c2
velf = df*pos + dg*vel

# compute the partials (Battin)
C = gm*(x^2)*(gm*(x^3)*(3c5 - c4) - dt*c2)
C = -3(x^3)*c3 + gm*(x^3) - gm*dt*(x^2)*c2
delr = posf - pos
delv = velf - vel

dxdx = f*Kepler.I3 + (rf/gm)*delv*transpose(delv) + (1/r0^3)*(r0*(1-f)*posf*transpose(pos) + C*velf*transpose(pos))
dxdx .- dxdx_auto
# -dxdx .- dxdx_auto