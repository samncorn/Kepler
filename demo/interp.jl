# # interpolate the stumpff function
using Pkg; Pkg.activate()
using GLMakie
using Kepler
# using DataInterpolations
using FastChebInterp

# function stumpff_interpolate(z)
    
# end

# function make_stumpff_table(::T) where {T}

# end
# zn = 10 .^ (-3:0.1:1)
# z_nodes  = -10:1e-3:10
# z_nodes = vcat(
#     1e-3:1e-3:9e-3,
#     1e-2:1e-3:9e-2,
#     1e-1:1e-2:9e-1,
#     1e-0:1e-1:9e-0,
#     1e+1:1e-1:9e+1,
# )

# m = 2^9
# z_nodes  = (2pi/m:2pi/m:sqrt(20)) .^ 2
# z_nodes  = vcat(-reverse(z_nodes), 0.0, z_nodes)

z_min = 0.0
z_max = 100.0
n_points = 2000
stumpff_c2(z) = Kepler.stumpff5(z)[3]
stumpff_c3(z) = Kepler.stumpff5(z)[4]
x = chebpoints(n_points, z_min, z_max)
c2_interp = chebinterp(stumpff_c2.(x), z_min, z_max)


# dc2_nodes = -c4_nodes .- (c3_nodes .- 4*c4_nodes) ./ 2
# c2_interp = CubicHermiteSpline(dc2_nodes, c2_nodes, z_nodes)
# # c2_interp = 

# dc3_nodes = -c5_nodes .- (c4_nodes .- 5*c5_nodes) ./ 2
# c3_interp = CubicHermiteSpline(dc3_nodes, c3_nodes, z_nodes)


# c2, c3
z  = 1e-2:1e-2:z_max
c_vals = Kepler.stumpff5.(z)
c0_vals = [x[1] for x in c_vals]
c1_vals = [x[2] for x in c_vals]
c2_vals = [x[3] for x in c_vals]
c3_vals = [x[4] for x in c_vals]

c22 = c2_interp.(z)
# c32 = c3_interp.(z)

c2_err = c2_vals .- c22
# c3_err = c3 .- c32

f  = Figure()
ax = Axis(f[1, 1])
# lines!(ax, z, c22)
lines!(ax, z, c2_err)
# lines!(ax, z, c3_err)
# lines!(ax, z, c0)
# # lines!(ax, z, c3)
f