module Kepler

using Roots
using StaticArrays
using LinearAlgebra
# using Rotations
# using Enzyme
using Printf

""" Everything needed to paramterize a keplerian orbit wrapped up nice
"""
struct Orbit{T, V}
    position::V
    velocity::V
    t::T
    gm::T
end

""" RADIANS """
struct CometaryOrbit{T}
    q::T
    e::T
    i::T
    Om::T
    w::T
    tp::T
    epoch::T
    gm::T
end

function vec_angle(v1, v2)
    n1 = norm(v1)
    n2 = norm(v2)
    return 2atan(norm(n1*v2 - n2*v1), norm(n1*v2 + n2*v1))
end

const I3 = SMatrix{3, 3}((
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    ))

xrot(x) = SMatrix{3, 3}((
    1.0, 0.0, 0.0,
    0.0, cos(x), sin(x),
    0.0, -sin(x), cos(x)
    ))

zrot(x) = SMatrix{3, 3}((
    cos(x), sin(x), 0.0,
    -sin(x), cos(x), 0.0,
    0.0, 0.0, 1.0,
    ))

include("roots.jl")
include("stumpff.jl")
include("propagate.jl")
include("moid.jl")
include("Lambert.jl")

# dot(v1, v2) = sum(v1 .* v2)
# norm(v) = sqrt(dot(v, v))
# cross(v1::SVector{3, T}, v2::SVector{3, T}) where {T} = T(
#     v1[2]*v2[3] - v1[3]*v2[2],
#     v1[3]*v2[1] - v1[1]*v2[3],
#     v1[1]*v2[2] - v1[2]*v2[1],
#     )
# cross(v1, v2) = cross(SVector{3}(v1), SVector{3}(v2))



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

    # I = SVector{3}()
    # angles
    # i = acos(hvec[3]/h)
    i = atan(sqrt(hvec[1]^2 + hvec[2]^2), hvec[3])

    # Om = acos(node[1])
    Om = atan(sqrt(node[2]^2 + node[3]^2), node[1])
    if node[2] < 0
        Om = 2pi - Om
    end

    # w = acos(dot(node, evec)/e)
    # w = atan()
    w = vec_angle(node, evec)
    if evec[3] < 0
        w = 2pi - w
    end

    # nu = acos(dot(evec, pos)/(e*r))
    nu = vec_angle(evec, pos)
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
    vq  = sqrt(gm*(1 + e)/q)
    pos = SVector{3}(q, 0.0, 0.0)
    vel = SVector{3}(0.0, vq, 0.0)

    # rotate 
    R  = zrot(Om)*xrot(i)*zrot(w)
    
    pos = R*pos
    vel = R*vel
    # propagate to epoch
    return propagate(pos, vel, epoch - tp, gm)
end

# function periapse(orbit)

# end

end # module