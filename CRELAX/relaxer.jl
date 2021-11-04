using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
include("MoS2_0_Parameters.jl") # replace with appropriate parameters file
using Optim
using LinearAlgebra
using NPZ
using ArgParse
# using PyPlot, PyCall

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

function main()
    ## Parameters:
    parsed_args = parse_commandline()
    dir = abspath(parsed_args["out"])
    println("Output directory: $(dir)")
    theta = deg2rad(parsed_args["deg"])
    println("Using twist angle $(theta) rad or $(rad2deg(theta)) deg")

    N = parsed_args["N"]
    println("Using $(N) x $(N) grid")
    hN = div(N,2)
    blg = Bilayer(l, theta, K, G)
    hull = Hull(blg, N);

    ## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
    f(u) = Energy(u, hull)
    g!(storage, u) = Gradient!(storage, u, hull)
    @time results = optimize(f, g!, zeros(2, N^2), LBFGS())
    # u is the relaxation of L1 wrt pristine L2, so if adding on to a configuration in L2
    # the sign flips and it becomes b <- b - u
    u = -2*reshape(Optim.minimizer(results), (2, N^2))
    bprime = transpose(u) + transpose(hull.G) # output
    npzwrite(dir*"/b_cart.npz", transpose(hull.G))
    npzwrite(dir*"/u_cart.npz", transpose(u))
    npzwrite(dir*"/bprime_cart.npz", bprime)
    println("Output written to $(parsed_args["out"])")
    if false
        GSFE_range = GSFE(blg.E*[2/3 1/3 .46385 0 ;2/3 1/3 .46385 0], blg)
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
        C = complete(reshape(GSFE(hull.G + hull.bl.iR * (Reflection(u, hull) - u), hull.bl), (N,N)));
        close("all")
        fig = figure(1, figsize = (8, 10))
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        clf()
        subplot(2,1,1)
        pcolormesh(X, Y, C0, cmap=MyColorMap, shading="gouraud")
        clim(GSFE_range[[1,end]])
        ax = gca()
        ax.set_aspect(1)
        subplot(2,1,2)
        pcolormesh(X, Y, C, cmap=MyColorMap, shading="gouraud")
        clim(GSFE_range[[1,end]])
        ax = gca()
        ax.set_aspect(1)
        savefig(dir*"/rgsfe.png")
        figure(1)
        println(max(u[1,:]...))
        println("Plotting done")
    end
end
main()

