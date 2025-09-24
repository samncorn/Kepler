using Kepler
using StaticArrays
using ProgressMeter

# sample a range of orbits
a = 0.01:0.01:100.0  # semimajor axis (au)
e = 0:0.01:10        # eccentricity
f = 0:(2pi/10):2pi # true anomaly (rad)
T = 10 .^ (-2:0.5:1)     # number of revolutions ()

println(log10(length(Iterators.product(a, e, f, T))))
@showprogress for (ai, ei, fi, Ti) in Iterators.product(a, e, f, T)
    # build initial state

    # compute period, propagation timespan

    # try to propagate, if fails, dump the initial conditions to run through a debug build
end