## Example of Graphene optimization + pattern plot in configuration space

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
#include("MoS2_180_Parameters.jl")
include("MoS2_0_Parameters.jl")
include("FourierTools.jl")
using PyPlot
using PyCall
using Interpolations
using Optim
using Polynomials
using LinearAlgebra # for dot
using Printf # simple macros for printing formatted output
using DelimitedFiles # writedlm, for file output

## Parameters:
#theta_list = linspace(0.3,2,11)
#theta_list = linspace(0.9, 1.3, 21)

theta_list = [0.5]
#theta_list = range(0.3,stop=2,length=15)

# Rotate the effective basis-vectors for the supercell
# Just rotates u_x and u_y into one another
# convention is CCW rotation for positive angles
basis_rot = deg2rad(-30.0);

max_p = 12;
u_x_coeffs = zeros(length(theta_list),2*max_p+1,2*max_p+1);
u_y_coeffs = zeros(length(theta_list),2*max_p+1,2*max_p+1);
u_z_coeffs = zeros(length(theta_list),2*max_p+1,2*max_p+1);
#x_cos_coeffs = zeros(length(theta_list),max_p,max_p);
#x_sin_coeffs = zeros(length(theta_list),max_p,max_p);
#y_cos_coeffs = zeros(length(theta_list),max_p,max_p);
#y_sin_coeffs = zeros(length(theta_list),max_p,max_p);

#coeffs_mag = zeros(0.0)
coeffs_xr = zeros(0)
coeffs_yr = zeros(0)
coeffs_xi = zeros(0)
coeffs_yi = zeros(0)
#coeffs_z = zeros(0.0)

coeffs_x = zeros(0)
coeffs_y = zeros(0)

k_list = zeros(0)

for theta_idx = 1:length(theta_list)
    θ = deg2rad(theta_list[theta_idx])

    N = 32
    hN = div(N,2)
    blg = Bilayer(l, θ, K, G)
    hull = Hull(blg, N);


    ## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
    f(u) = Energy(u, hull)
    g!(storage, u) = Gradient!(storage, u, hull)
    @time results = optimize(f, g!, zeros(2, N^2), LBFGS())

    u = reshape(Optim.minimizer(results), (2, N^2))
    ∇u = ∇(u, hull)


    GSFE_range = GSFE(blg.E*[2/3 1/2 0 ;2/3 1/2 0], blg)
    mid = (GSFE_range[2]-GSFE_range[1])/(GSFE_range[3]-GSFE_range[1])
    cscale(c) = c < mid ? .5*c/mid : .5 *(1. + (c-mid)/(1-mid))
    MyColorMap = ColorMap( ColorMap("inferno")(cscale.(range(0,stop=1.,length=512))))

    function complete(A::Array{Float64,2})
        A = hcat(A, A[:,1])
        A = vcat(A, A[1,:]')
    end
    function complete(A::Array{Float64,2}, hoffset::Float64, voffset::Float64)
        A = hcat(A, A[:,1].+hoffset)
        A = vcat(A, A[1,:]'.+voffset)
    end


    # Gather and shape output
    X = complete(reshape(hull.G[1,:], (N,N)), dot(blg.E[1,:],[0;1]), dot(blg.E[1,:],[1;0]))
    Y = complete(reshape(hull.G[2,:], (N,N)), dot(blg.E[2,:],[0;1]), dot(blg.E[2,:],[1;0]))
    C0 = complete(reshape(GSFE(hull.G, hull.bl), (N,N)))
    C = complete(reshape(GSFE(hull.G + hull.bl.iR * (Reflection(u, hull) - u), hull.bl), (N,N)))
    #Z = reshape(GSFE_height(hull.G + hull.bl.iR * (Reflection(u, hull) - u), hull.bl), (N,N))

    #U_x_temp = (reshape(u[1,:], (N,N)))
    #U_y_temp = (reshape(u[2,:], (N,N)))
    u_fixed = (hull.bl.iR * (Reflection(u, hull) - u))
    U_x_temp = reshape( u_fixed[1,:] , (N,N) )
    U_y_temp = reshape( u_fixed[2,:] , (N,N) )

    #U_z_temp = Z
    config_x_arr = reshape(hull.G[1,:], (N,N))
    config_y_arr = reshape(hull.G[2,:], (N,N))

    R_t2 = [cos(θ/2.0) -sin(θ/2.0); sin(θ/2.0) cos(θ/2.0)]

    # config space to realspace converter
    cf2r = inv(I-R_t2*R_t2)*R_t2

    sc_x_arr = cf2r[1,1]*config_x_arr + cf2r[1,2]*config_y_arr
    sc_y_arr = cf2r[2,1]*config_x_arr + cf2r[2,2]*config_y_arr

    #=
    println(config_x_arr[end,1])
    println(config_y_arr[end,1])
    println(config_x_arr[1,end])
    println(config_y_arr[1,end])

    println(cf2r)

    println(sc_x_arr[end,1])
    println(sc_y_arr[end,1])
    println(sc_x_arr[1,end])
    println(sc_y_arr[1,end])
    =#

    sc_E = cf2r*blg.E

    #println(size(u))
    #println(size(U_x_temp))

    # put relaxation into unit-cell basis

    #(3, 2, 2) 0.000000 + 0.000593 im

    U_X = cos(basis_rot)*U_x_temp - sin(basis_rot)*U_y_temp
    U_Y = sin(basis_rot)*U_x_temp + cos(basis_rot)*U_y_temp
    #U_Z = U_z_temp

    k_list = fourierGetSamplingFull(sc_E, max_p)

    global u_x_coeffs = fourierDecompose(U_X/2.0,sc_x_arr,sc_y_arr,k_list)
    global u_y_coeffs = fourierDecompose(U_Y/2.0,sc_x_arr,sc_y_arr,k_list)
    #u_z_coeffs = fourierDecompose((U_Z-3.35)/2.0,sc_x_arr,sc_y_arr,k_list)
    #mags = sqrt.(abs(u_x_coeffs).^2 + abs(u_y_coeffs).^2)

    #=
    println("u_x")
    show(STDOUT, "text/plain", imag(u_x_coeffs))
    println("\nu_y")
    show(STDOUT, "text/plain", imag(u_y_coeffs))
    println("\nu_z")
    show(STDOUT, "text/plain", real(u_z_coeffs))
    =#

    #=
    if (theta_idx == 1)
        coeffs_mag = mags;
    else
        coeffs_mag = hcat(coeffs_mag,mags)
    end
    =#
    if (theta_idx == 1)
        #println("u_x :", u_x_coeffs)
        #println("u_y :", u_y_coeffs)

        global coeffs_xr = real(u_x_coeffs)
        global coeffs_yr = real(u_y_coeffs)
        global coeffs_xi = imag(u_x_coeffs)
        global coeffs_yi = imag(u_y_coeffs)
        #coeffs_z = real(u_z_coeffs)
    else
        global coeffs_x = hcat(coeffs_x,imag(u_x_coeffs))
        global coeffs_y = hcat(coeffs_y,imag(u_y_coeffs))
        #coeffs_z = hcat(coeffs_z,real(u_z_coeffs))
    end
end

if length(theta_list) == 1
    for k_idx = 1:size(k_list,2)
        #println(  @sprintf("u_x [%d, %d]: %.10f", k_list[3,k_idx], k_list[4,k_idx], coeffs_x[k_idx]) )
    end
    for k_idx = 1:size(k_list,2)
        #println(  @sprintf("u_y [%d, %d]: %.10f", k_list[3,k_idx], k_list[4,k_idx], coeffs_y[k_idx]) )
    end
    for k_idx = 1:size(k_list,2)
        #println(  @sprintf("u_z [%d, %d]: %.10f", k_list[3,k_idx], k_list[4,k_idx], coeffs_z[k_idx]) )
    end

    # Write data to file
    # #=
    writedlm("data/thetas.txt", theta_list, ", ")
    writedlm("data/coeffs_xr.txt", coeffs_xr, ", ")
    writedlm("data/coeffs_yr.txt", coeffs_yr, ", ")
    writedlm("data/coeffs_xi.txt", coeffs_xi, ", ")
    writedlm("data/coeffs_yi.txt", coeffs_yi, ", ")
    #writedlm("data/0p1_2p0_k15/coeffs_z.txt", coeffs_z, ", ")
    #writedlm("data/0p1_2p0_k15/elastic_energies.txt", elastic_E, ", ")

    # =#


else

    #close("all")
    fig = figure(1, figsize = (8, 10))

    subplot2grid((3, 1), (0, 0))
    for k_idx = 1:size(k_list,2)
        plot(theta_list,coeffs_x[k_idx,:],"-x")
        #semilogy(theta_list,coeffs_mag[k_idx,:])
        hold(true)
        label = @sprintf("[%d, %d]",k_list[3,k_idx],k_list[4,k_idx])
        text(theta_list[end]-.03,coeffs_x[k_idx,end]*1.1,label)

        p = polyfit(theta_list-theta_list[end],coeffs_x[k_idx,:],2)
        plot(theta_list,p(theta_list-theta_list[end]),"--k")
        println(  @sprintf("u_x [%d, %d]: %s", k_list[3,k_idx], k_list[4,k_idx], coeffs(p)) )

    end

    xlabel("Twist (deg)")
    ylabel("U_X (Ang)")

    subplot2grid((3, 1), (1, 0))
    for k_idx = 1:size(k_list,2)
        plot(theta_list,coeffs_y[k_idx,:],"-x")
        #semilogy(theta_list,coeffs_mag[k_idx,:])
        hold(true)
        label = @sprintf("[%d, %d]",k_list[3,k_idx],k_list[4,k_idx])
        text(theta_list[end]-.03,coeffs_y[k_idx,end]*1.1,label)

        p = polyfit(theta_list-theta_list[end],coeffs_y[k_idx,:],2)
        plot(theta_list,p(theta_list-theta_list[end]),"--k")
        println(  @sprintf("u_y [%d, %d]: %s", k_list[3,k_idx], k_list[4,k_idx], coeffs(p)) )

    end



    xlabel("Twist (deg)")
    ylabel("U_Y (Ang)")

    subplot2grid((3, 1), (2, 0))
    for k_idx = 1:size(k_list,2)
        #plot(theta_list,coeffs_z[k_idx,:],"-x")
        #semilogy(theta_list,coeffs_mag[k_idx,:])
        #hold(true)
        #label = @sprintf("[%d, %d]",k_list[3,k_idx],k_list[4,k_idx])
        #text(theta_list[end]-.03,coeffs_z[k_idx,end]*1.1,label)

        #p = polyfit(theta_list-theta_list[end],coeffs_z[k_idx,:],2)
        #plot(theta_list,p(theta_list-theta_list[end]),"--k")
        #println(  @sprintf("u_z [%d, %d]: %s", k_list[3,k_idx], k_list[4,k_idx], coeffs(p)) )

    end

    #xlabel("Twist (deg)")
    #ylabel("U_Z (Ang)")

    #=
    for k_idx = 1:size(k_list,2)
        #plot(theta_list,coeffs_mag[k_idx,:])
        semilogy(theta_list,coeffs_mag[k_idx,:])
        hold(true)
        label = @sprintf("[%d, %d]",k_list[3,k_idx],k_list[4,k_idx])
        text(theta_list[end]-.03,coeffs_mag[k_idx,end]*1.1,label)
    end
    =#

    hold(false)
end
