## Example of Graphene optimization + pattern plot in configuration space
using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)

include("Bilayers.jl")
include("gsfe_func.jl")

using Interpolations
using Optim
using LinearAlgebra
using ArgParse
using NPZ

on_cluster = false #! make true when on cluster

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--deg", "-d"
            help = "twist angle in degrees"
            required = true
            arg_type = Float64
        "-N", "-n"
            help = "grid size (N if grid is NxN)"
            arg_type = Int
            required = true
        "--out", "-o"
            help = "output directory path (must end in slash)"
            arg_type = String
            default = "."
    end

    return parse_args(s)
end

## Parameters:
parsed_args = parse_commandline()
dir = string(abspath(parsed_args["out"]), "/")
println("Output directory: $(dir)")
θ = deg2rad(parsed_args["deg"])

N = parsed_args["N"]
println("Using $(N) x $(N) grid")
hN = div(N,2)
blg = Bilayer(l, θ, K1, G1, K2, G2)
hull = Hull(blg, N, [0., 0.]);


## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
f(u) = Energy(u, hull)
g!(storage, u) = Gradient!(storage, u, hull)
@time results = optimize(f, g!, zeros(4, N^2), LBFGS(),
    Optim.Options(allow_f_increases=true, show_trace = true, 
                  store_trace = true, show_every = 100))
println(results)

u = reshape(Optim.minimizer(results), (4, N^2))
u1 = u[1:2, :]
u2 = u[3:4, :]
rel_u = u1 + u2 # relative relaxation between two layers is twice as strong
GSFE_range = GSFE(blg.E*[2/3 1/3 .46385 0 ;2/3 1/3 .46385 0], blg)
GSFE_range = sort(GSFE_range)

rel_u_T = transpose(rel_u)
b_cart = transpose(hull.G)
npzwrite(dir*"b_cart.npz", b_cart)
npzwrite(dir*"u_cart.npz", rel_u_T)
npzwrite(dir*"bprime_cart.npz", rel_u_T + b_cart)
println("Output written to $(parsed_args["out"])")

if !on_cluster
    using PyPlot
    using PyCall
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
    C = complete(reshape(GSFE(hull.G + hull.bl.iR * (u1 + u2), hull.bl), (N,N)));

    # convert to supercell positions
    r = inv(blg.R - blg.iR) * [X[:], Y[:]]
    x = r[1]
    y = r[2]
    ux_l1 = complete(reshape(u[1,:], (N,N)))
    uy_l1 = complete(reshape(u[2,:], (N,N)))
    ux_l2 = complete(reshape(u[3,:], (N,N)))
    uy_l2 = complete(reshape(u[4,:], (N,N)))


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
    savefig(string( dir, "GSFE_" , rad2deg(θ), ".png"))
    figure(1)

    figure(2)
    quiver(reshape(x,(N+1,N+1)), reshape(y, (N+1,N+1)), ux_l1, uy_l1)
    ax = gca()
    ax.set_aspect(1)
    savefig(string( dir, "RelField1_" , rad2deg(θ), ".png"))

    figure(3)
    quiver(reshape(x,(N+1,N+1)), reshape(y, (N+1,N+1)), ux_l2, uy_l2)
    ax = gca()
    ax.set_aspect(1)
    savefig(string( dir, "RelField2_" , rad2deg(θ), ".png"))


    u1x = reshape(ux_l1, (N+1, N+1))
    u1y = reshape(uy_l1, (N+1, N+1))
    u2x = reshape(ux_l2, (N+1, N+1))
    u2y = reshape(uy_l2, (N+1, N+1))
    println(string("Layer 1 max ux: ", max(u1x...)))
    println(string("Layer 2 max ux: ", max(u2x...)))
end