""" I believe this to be the proper declension
"""
abstract type Osculator{T} end

#--- Cartesian State ---
struct Cartesian{T, V} <: Osculator{T}
    position::V
    velocity::V
    epoch::T
    gm::T
end

struct Cometary{T} <: Osculator{T}
    q::T
    e::T
    i::T
    Om::T
    w::T
    tp::T
    epoch::T
    gm::T
end

struct Keplerian{T} <: Osculator{T}
    a::T
    e::T
    i::T
    Om::T
    w::T
    M::T
    epoch::T
    gm::T
end

function Cartesian(elements::Cometary)
    q = elements.q
    e = elements.e
    i = elements.i
    Om = elements.Om
    w = elements.w
    tp = elements.tp
    epoch = elements.epoch
    gm = elements.gm
    # compute periapse pos, vel in perifocal frame
    vq  = sqrt(gm*(1 + e)/q)
    pos = SVector{3}(q, 0.0, 0.0)
    vel = SVector{3}(0.0, vq, 0.0)

    # rotate 
    R  = zrot(-Om)*xrot(-i)*zrot(-w)
    
    pos = R*pos
    vel = R*vel
    # propagate to epoch
    pos, vel = propagate(pos, vel, epoch - tp, gm)
    return Cartesian(pos, vel, epoch, gm)
end

Cartesian(elements::Keplerian) = Cartesian(Cometary(elements))

# """ Get the variational equations  
# """
# function variations(state::Cartesian, t)

# end

#--- Keplerian Elements ---

Keplerian(state::Cartesian) = Keplerian(Cometary(state))

function Keplerian(elements::Cometary)
    a = elements.q/(1.0 - elements.e)
    n = sqrt(gm/elements.a)
    M = elements.tp*n
    return Keplerian(a, elements.e, elements.i, elements.Om, elements.w, M, elements.epoch, elements.gm)
end

# need to iron out interface
# function variations(P, elements::Keplerian)
#     da  = 
#     de  =
#     di  = 
#     dOm =
#     dw  = 
#     dM  =
#     return da, de, di, dOm, dw, dM
# end

#--- Cometary Elements ---

function Cometary(elements::Keplerian)
    q  = elements.a * (1.0 - elements.e) 
    n  = sqrt(gm/elements.a^3)
    tp = elements.M/n
    return Cometary(q, elements.e, elements.i, elements.Om, elements.w, tp, elements.epoch, elements.gm)
end

function Cometary(state::Cartesian)
    pos   = state.position
    vel   = state.velocity
    gm    = state.gm
    epoch = state.epoch

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

    # return q, e, i, Om, w, epoch - dt
    return Cometary(q, e, i, Om, w, dt, epoch, gm)
end

# TODO
#--- Equinoctial Elements ---
# struct ModifiedEquinoctialState{T}

# end

