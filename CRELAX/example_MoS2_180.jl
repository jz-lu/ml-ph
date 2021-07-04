using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
include("GrapheneParameters.jl")
# using PyPlot
using PyCall
using Interpolations
using Optim
using LinearAlgebra

## Parameters:
θ = deg2rad(1)

N = 50
hN = div(N,2)
blg = Bilayer(l, θ, K, G)
hull = Hull(blg, N);


## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
f(u) = Energy(u, hull)
g!(storage, u) = Gradient!(storage, u, hull)
@time results = optimize(f, g!, zeros(2, N^2), LBFGS())
u = reshape(Optim.minimizer(results), (2, N^2))
print(u)
# open(string(".", round(rad2deg(θ), digits=3), "u.txt"), "w") do f
#     write(f, "x,y,u1x,u1y,u2x,u2y \n")
#     for i in 1:(N+1)^2
#         xtemp = x[i]
#         ytemp = y[i]
#         r1 = u1x[i]
#         r2 = u1y[i]
#         r3 = u2x[i]
#         r4 = u2y[i]
#         write(f, "$xtemp,$ytemp,$r1,$r2,$r3,$r4\n")
#     end
# end

# ∇u = ∇(u, hull)

# GSFE_range = GSFE(blg.E*[2/3 1/3 .46385 0 ;2/3 1/3 .46385 0], blg)
# mid1 = (GSFE_range[2]-GSFE_range[1])/(GSFE_range[4]-GSFE_range[1])
# mid2 = (GSFE_range[3]-GSFE_range[1])/(GSFE_range[4]-GSFE_range[1])
# cscale(c) = c < mid1 ? c/(3*mid1) : (c < mid2 ? ( 1. + (c-mid1)/(mid2-mid1))/3. : (2. + (c-mid2)/(1-mid2))/3.)
# MyColorMap = ColorMap( ColorMap("Spectral_r")(cscale.(range(0,stop=1.,length=512))))


# function complete(A::Array{Float64,2})
#     A = hcat(A, A[:,1])
#     A = vcat(A, A[1,:]')
# end
# function complete(A::Array{Float64,2}, hoffset::Float64, voffset::Float64)
#     A = hcat(A, A[:,1].+hoffset)
#     A = vcat(A, A[1,:]'.+voffset)
# end


# X = complete(reshape(hull.G[1,:], (N,N)), dot(blg.E[1,:],[0;1]), dot(blg.E[1,:],[1;0]))
# Y = complete(reshape(hull.G[2,:], (N,N)), dot(blg.E[2,:],[0;1]), dot(blg.E[2,:],[1;0]))
# C0 = complete(reshape(GSFE(hull.G, hull.bl), (N,N)));
# C = complete(reshape(GSFE(hull.G + hull.bl.iR * (Reflection(u, hull) - u), hull.bl), (N,N)));


# close("all")
# fig = figure(1, figsize = (8, 10))
# fig.subplots_adjust(hspace=0.1, wspace=0.1)
# clf()
# subplot(2,1,1)
# pcolormesh(X, Y, C0, cmap=MyColorMap, shading="gouraud")
# clim(GSFE_range[[1,end]])
# ax = gca()
# ax.set_aspect(1)

# subplot(2,1,2)
# pcolormesh(X, Y, C, cmap=MyColorMap, shading="gouraud")
# clim(GSFE_range[[1,end]])
# ax = gca()
# ax.set_aspect(1)
# figure(1)


# print(max(u[1,:]...))
