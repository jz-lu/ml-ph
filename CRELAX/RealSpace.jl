## Given displacement discretized on a bilayer Hull, 
# provide functions to compute relaxed and unrelaxed 
# real space positions as well as resulting diffraction patterns.

include("GrapheneBilayers.jl")
using Interpolations
mutable struct Configuration
    γ::Array{Float64,2}

    function Configuration(γ1::Array{Float64,1}, γ2::Array{Float64,1}, hull::Hull)
        γ = hcat(   hull.bl.tE[1]*mod.(inv(hull.bl.tE[1])*γ1, 1.), 
                    hull.bl.tE[2]*mod.(inv(hull.bl.tE[2])*γ2, 1.))
        new(γ)
    end
end

mutable struct Displacement
    itp::Array{Interpolations.Extrapolation,2}     # Array of interpolants for the modulation functions

    function Displacement(u::Array{Float64,2}, hull::Hull)
        N = hull.N
        itp = Array{Interpolations.Extrapolation, 2}(2,2)

        function periodize(A::Array{Float64,2})
            A = hcat(A, A[:,1])
            A = vcat(A, A[1,:]')
        end
        U = periodize(reshape(u[1,:,:], (N,N)))
        itp[1,1] = extrapolate(interpolate(U, BSpline(Cubic(Line())), OnGrid()), Periodic())
        U = periodize(reshape(u[2,:,:], (N,N)))
        itp[2,1] = extrapolate(interpolate(U, BSpline(Cubic(Line())), OnGrid()), Periodic())

        # Reflection u -> S u(S gamma)
        u2 = permutedims(reshape(u, (2,N,N)), [1,3,2])
        u2[2,:,:] *= -1.
        U = periodize(reshape(u2[1,:,:], (N,N)))
        itp[1,2] = extrapolate(interpolate(U, BSpline(Cubic(Line())), OnGrid()), Periodic())
        U = periodize(reshape(u2[2,:,:], (N,N)))
        itp[2,2] = extrapolate(interpolate(U, BSpline(Cubic(Line())), OnGrid()), Periodic())

        new(itp)
    end
end


# Generate the unrelaxed positions of atoms in a given radius for a given configuration

function create_positions(R::Float64, ω::Configuration, hull::Hull)
    n = floor(Int64, R/norm(inv(hull.bl.E)))+1
    I = (-n:n) .* ones(Int64, (1,2*n+1))
    J = ones(Int64, (2*n+1,1)) .* (-n:n)'

    X1 = ω.γ[:,1] .+ hull.bl.tE[1] * [I[:]'; J[:]'];
    X2 = ω.γ[:,2] .+ hull.bl.tE[2] * [I[:]'; J[:]'];

    X1 = X1[:, hypot.(X1[1,:],X1[2,:]) .< R]
    X2 = X2[:, hypot.(X2[1,:],X2[2,:]) .< R]
    return [[X1]; [X2]]
end

# Relax the position of atoms in both layers according to the provided displacement.
# Positions is an array of arrays of unrelaxed positions indexed by the layer index (1,2), 
#    where the first dimension corresponds to directions x and y,
#          the second dimension to point index.
function displace!(Positions::Array{Array{Float64, 2}, 1}, ϕ::Displacement, ω::Configuration, hull::Hull)
    # Layer 1
    CartPos = hull.N .* (inv(hull.bl.tE[2]) * (ω.γ[:,2] .- Positions[1])).+1
    for i=1:2, ind=1:size(Positions[1], 2)
        Positions[1][i,ind] += ϕ.itp[i,1][CartPos[1,ind], CartPos[2,ind]]
    end
    # Layer 2
    CartPos = hull.N .* (inv(hull.bl.tE[1]) * (ω.γ[:,1] .- Positions[2])) .+ 1
    for i=1:2, ind=1:size(Positions[2], 2)
        Positions[2][i,ind] += ϕ.itp[i,2][CartPos[1,ind], CartPos[2,ind]]
    end
    nothing
end





function fourier_zoom(Positions::Array{Array{Float64,2},1}, 
                        sigma::Float64, 
                        K::Array{Float64,1}, 
                        L::Float64, 
                        N::Int64 = 64., 
                        m::Float64 = 1.)
    grid = (-N:N)

    filtereddata = zeros(Complex128, 2*N+1, 2*N+1)
    X = [0.;0.]
    for j=1:length(Positions)
        for index = 1:size(Positions[j],2)
            X[1] = Positions[j][1,index]/L;
            X[2] = Positions[j][2,index]/L;
            ϕ = exp(-2im*pi*L*dot(K,X) -.5*(L/sigma)^2*dot(X,X))
            for k = max(1, N+1 + floor(Int64, X[1]) - 6):min(2*N, N+1 + ceil(Int64, X[1]) + 6),
                l = max(1, N+1 + floor(Int64, X[2]) - 6):min(2*N, N+1 + ceil(Int64, X[2]) + 6)
                filtereddata[k,l] += ϕ * exp(-.5*((grid[k] - X[1])^2  + (grid[l] - X[2])^2))
            end
        end
    end

    filter = zeros(2*N+1)
    for k=max(1, N-6):min(2*N+1, N+6)
        filter[k] = exp(-.5*grid[k]^2)
    end
    intensity = abs2.(fftshift(fft(filter)))

    S = sum([sum(exp.( (-.5/sigma^2) .* sum(Pos.^2, 1))) for Pos in Positions])
    farfieldintensity = abs2.(fftshift(fft(filtereddata))) ./ intensity ./ intensity' ./ S


    zoom = 1+floor(Int64, 2*N/3):ceil(Int64, 4*N/3)
    Kx = K[1] .+ (m/(L*(2*N+1))) .* grid[zoom]' .* ones(zoom)
    Ky = K[2] .+ (m/(L*(2*N+1))) .* ones(zoom') .* grid[zoom]
    I = farfieldintensity[zoom,zoom]'

    return Kx, Ky, I
end
