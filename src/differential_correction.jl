# don't over complicate it. There are a finit number of obervable types, so we define a container that holds (potentially empty) containers of each
# Then we can define
# struct ObervablesSet{T}
#     optical::Vector{OpticalObservable{T}}
#     # radar::Vector{RadarObservable{T}}
# end

# function compute_residuals(obs::ObervablesSet{T}, x) where {T}
#     return Iterators.flatten(Iterators.map(((o, x),) -> compute_residuals(o, x), ))
# end


""" orbit fitting on optical measurements.
"""
function fit_orbit(compute_residuals <: Function, pos0, vel0; kwargs)

end

""" Orbit fitting with RANSAC-style outlier detection. Uses Herget's method for IOD, but still requires an initial guess to get initial ranges.
"""
function RANSAC(compute_residuals, pos0, vel0; n_samples = 1, kwargs)

end

# how to incoporate other meas types? i.e. optical + radar + ddor + occultation
# outline (OD loop)
#   x -> traj
#   for each meas
#       