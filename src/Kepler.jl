module Kepler

using Roots
using StaticArrays
using LinearAlgebra
using Enzyme
using Printf

# Write your package code here.
include("roots.jl")
include("stumpff.jl")
include("propagate.jl")
include("Lambert.jl")

xrot(x) = SVector{3, 3}((
    1.0, 0.0, 0.0,
    0.0, cos(x), sin(x),
    0.0, -sin(x), cos(x)
    ))

zrot(x) = SVector{3, 3}((
    cos(x), sin(x), 0.0,
    -sin(x), cos(x), 0.0,
    0.0, 0.0, 1.0,
    ))

"""returns radians"""
function cometary(pos, vel, epoch, gm)
    # shape params
    hvec = cross(pos, vel)
    node = normalize(cross(SVector{3}(0.0, 0.0, 1.0), hvec))
    v2   = dot(vel, vel)
    r    = norm(pos)
    evec = ((v2 - gm/r)*pos - dot(pos, vel)*vel)/gm
    e    = norm(evec)
    h2   = dot(hvec, hvec)
    h    = sqrt(h2)
    q    = h2/gm/(1 + e)
    l    = (1-e)/(1+e)

    # angles
    i = acos(hvec[3]/h)

    Om = acos(node[1])
    if node[2] < 0
        Om = 2pi - Om
    end

    w = acos(dot(node, evec)/e)
    if evec[3] < 0
        w = 2pi - w
    end

    nu = acos(dot(evec, pos)/(e*r))
    if dot(pos, vel) < 0
        nu = 2pi - nu
    end
    tau = tan(nu/2)
    # should make a better parablic case computation to smooth numeric problems
    # series given in Fukushima (1999)
    y = if l > 0
        2atan(sqrt(l)*tau)/sqrt(l)
    elseif l < 0
        2atanh(sqrt(-l)*tau)/sqrt(-l)
    else
        2tau
    end

    # times
    v = sqrt(gm*(1 + e)/(q^3))
    k1 = 1.0
    # k2 = dr0/(v*(q^2))
    k2 = 0.0
    k3 = gm/((v^2)*(q^3))
    L  = universal_kepler(y, l, k1, k2, k3)
    dt = L/v

    return q, e, i, Om, w, epoch - dt
end

""" assumes radians
"""
function cartesian(q, e, i, Om, w, tp, epoch, gm)
    # compute periapse pos, vel in perifocal frame
    pos = SVector{3}(q, 0.0, 0.0)
    vel = SVector{3}(0.0, sqrt(gm*(1+e)/q), 0.0)

    # rotate 
    R  = zrot(-Om)*xrot(-i)*zrot(-w)
    
    pos = R*pos
    vel = R*vel
    # propagate to epoch
    return propagate(pos, vel, epoch - tp, gm)
end

end # module