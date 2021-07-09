using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
include("MoS2_0_Parameters.jl") # replace with appropriate parameters file
using Optim
using LinearAlgebra
using NPZ
using ArgParse

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
    println("Using grid size $(N)")
    hN = div(N,2)
    blg = Bilayer(l, theta, K, G)
    hull = Hull(blg, N);

    ## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
    f(u) = Energy(u, hull)
    g!(storage, u) = Gradient!(storage, u, hull)
    @time results = optimize(f, g!, zeros(2, N^2), LBFGS())
    u = reshape(Optim.minimizer(results), (2, N^2))
    bprime = transpose(u) + transpose(hull.G) # output
    npzwrite(dir*"b_cart.npz", transpose(hull.G))
    npzwrite(dir*"u_cart.npz", transpose(u))
    npzwrite(dir*"bprime_cart.npz", bprime)
    println("Output written to $(parsed_args["out"])")
end
main()
