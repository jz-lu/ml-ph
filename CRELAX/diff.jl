## Example of Graphene optimization + pattern plot in configuration space
using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)

include("Bilayers.jl")
# include("MoS2_0_Parameters.jl")
include("MoS2_180_Parameters.jl")
#include("GrapheneParameters.jl")

using PyPlot
using PyCall
using Interpolations
using Optim
using LinearAlgebra

## Parameters:

θ = deg2rad(0.8)

N =54
hN = div(N,2)
blg = Bilayer(l, θ, K, G)
hull = Hull(blg, N, [0., 0.]);


## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
f(u) = Energy(u, hull)
g!(storage, u) = Gradient!(storage, u, hull)
@time results = optimize(f, g!, zeros(2, N^2), LBFGS(), Optim.Options(allow_f_increases=true))

u = reshape(Optim.minimizer(results), (2, N^2))

GSFE_range = GSFE(blg.E*[2/3 1/3 .46385 0 ;2/3 1/3 .46385 0], blg)
GSFE_range = sort(GSFE_range)
mid1 = (GSFE_range[2]-GSFE_range[1])/(GSFE_range[4]-GSFE_range[1])
mid2 = (GSFE_range[3]-GSFE_range[1])/(GSFE_range[4]-GSFE_range[1])
cscale(c) = c < mid1 ? c/(3*mid1) : (c < mid2 ? ( 1. + (c-mid1)/(mid2-mid1))/3. : (2. + (c-mid2)/(1-mid2))/3.)
MyColorMap = ColorMap( ColorMap("Spectral_r")(cscale.(range(0,stop=1.,length=512))))


function complete(A::Array{Float64,2})
    A = hcat(A, A[:,1])
    A = vcat(A, A[1,:]')
end
function complete(A::Array{Float64,2}, hoffset::Float64, voffset::Float64)
    A = hcat(A, A[:,1].+hoffset)
    A = vcat(A, A[1,:]'.+voffset)
end


X = complete(reshape(hull.G[1,:], (N,N)), dot(blg.E[1,:],[0;1]), dot(blg.E[1,:],[1;0]))
Y = complete(reshape(hull.G[2,:], (N,N)), dot(blg.E[2,:],[0;1]), dot(blg.E[2,:],[1;0]))
C0 = complete(reshape(GSFE(hull.G, hull.bl), (N,N)));
C = complete(reshape(GSFE(hull.G + hull.bl.iR * 2 * u, hull.bl), (N,N)));

# convert to supercell positions
r = inv(blg.R - blg.iR) * [X[:], Y[:]]
x = r[1]
y = r[2]
ux = complete(reshape(u[1,:], (N,N)))
uy = complete(reshape(u[2,:], (N,N)))


##
close("all")
fig = figure(1,figsize = (8, 10))
fig.subplots_adjust(hspace=0.1, wspace=0.1)
clf()
subplot(2,1,1)
pcolormesh(reshape(x,(N+1, N+1)), reshape(y, (N+1, N+1)), C0, cmap=MyColorMap, shading="gouraud")
clim(GSFE_range[[1,4]].-1e-3)
title("Unrelaxed")
ax = gca()
ax.set_aspect(1)
colorbar()

subplot(2,1,2)
pcolormesh(reshape(x,(N+1, N+1)), reshape(y, (N+1, N+1)), C, cmap=MyColorMap, shading="gouraud")
clim(GSFE_range[[1,4]].-1e-3)
ax = gca()
ax.set_aspect(1)
title("Relaxed")
colorbar()
savefig(string( "/Users/jonathanlu/Documents/twbl-phonon/Julia/my_figs/", "MoS2_0", "_" , rad2deg(θ), ".png"))
figure(1)

figure(2)
quiver(reshape(x,(N+1,N+1)), reshape(y, (N+1,N+1)), ux, uy )
ax = gca()
ax.set_aspect(1)
savefig(string( "/Users/jonathanlu/Documents/twbl-phonon/Julia/my_figs/", "MoS2_0_2", "_" , rad2deg(θ), ".png"))


u1x = reshape(ux, (N+1, N+1))
u1y = reshape(uy, (N+1, N+1))
print(max(u1x...))

