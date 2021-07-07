### Relaxation of bilayers using Hull framework.
## Common setup to all rotated homo-bilayers.

using LinearAlgebra
using NPZ

mutable struct Bilayer
    # Unrotated monolayer parameters
    l::Float64
    E::Array{Float64,2}
    invE::Array{Float64,2}
    P::Array{Float64,1}

    # Twist, rotation
    θ::Float64
    R::Array{Float64,2}
    iR::Array{Float64,2}

    # Twisted layers parameterization
    tE::Array{Array{Float64,2},1}
    tP::Array{Array{Float64,1},1}

    # Effective elastic tensor in hull space
    tC::Array{Float64,4}

    function Bilayer(l, θ, K, G)
        E = l * sqrt(3) * [1 1/2; 0 sqrt(3)/2]
        P = l * [1; 0]

        ## Symmetries and rotations:
        R  = [cos(θ/2) sin(θ/2); -sin(θ/2) cos(θ/2)]
        iR = [cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)]

        ## Lattice parameters for twisted layers:
        tE = [iR*E,R*E]
        tP = [iR*P,R*P]

        ## Intra-layer elastic tensor, with stress tensor linear ordering [E_11; E_12; E_21; E_22]
        C = [K+G 0 0 K-G; 0 G G 0; 0 G G 0; K-G 0 0 K+G]
        T = inv(E)  * (R - iR)
        T = [T[1,1] 0 T[1,2] 0; 0 T[1,1] 0 T[1,2]; T[2,1] 0 T[2,2] 0; 0 T[2,1] 0 T[2,2]] # Hull [0,1]^2 -> Real space gradient transformation
        tC =  T * C * T';
        tC = reshape(tC, (2,2,2,2))

        new(l, E, 2*pi*inv(E), P, θ, R, iR, tE, tP, tC)
    end
end

# Discretized Hull structure
mutable struct Hull
    bl::Bilayer
    N::Int
    G::Array{Float64,2}
    Gref::Array{Float64,2}
    Hess_Elastic::Array{Float64,4}
    plan
    iplan

    # Generate Hull discretization
    function Hull(bl::Bilayer, N)
        s = range(0,stop=1-1/N,length=N) .* ones(1,N)
        t = ones(N,1) .* range(0,stop=1-1/N,length=N)'
        G = bl.E*[s[:]'; t[:]']
        Gref = 2*pi*[s[:]'; t[:]']

        hN = div(N, 2)
        if iseven(N)  # even number of points
            K = 2*pi*[0:hN;1-hN:-1]
        else # odd number of points
            K = 2*pi*[0:hN; -hN:-1]
        end
        K1 = reshape(K, (N, 1))
        K2 = reshape(K, (1, N))

        A = zeros(Float64, (2, N, N))
        plan = plan_fft(A, (2,3) )
        A = zeros(Complex{Float64}, (2, N, N))
        iplan = plan_ifft!(A, (2,3))

        Hess = zeros(2,2,N,N)
        for a=1:2, c=1:2
            Hess[a,c,:,:] =  bl.tC[a,1,c,1] * (K1.*K1) .+ (bl.tC[a,2,c,1] .+ bl.tC[a,1,c,2]) * (K1.*K2) .+ bl.tC[a,2,c,2] * (K2.*K2)
        end

        new(bl, N, G, Gref, Hess, plan, iplan)
    end
end

# Implements reflection u -> S u(-S gamma)
function Reflection(u::Array{Float64, 2}, hull::Hull)
    s = reshape(u, (2, hull.N, hull.N))
    s = permutedims(s, [1,3,2])
    s[1,:,:] *= -1.
    return reshape(s, (2, hull.N^2))
end

# Implements GSFE energy for a displacement u discretized on the hull
function Misfit(u::Array{Float64, 2}, hull::Hull)
    return sum(GSFE(hull.G + hull.bl.iR * (Reflection(u, hull) - u), hull.bl))
    # return sum(GSFE(hull.G + hull.bl.iR * 2*u, hull.bl))
end

# Implements first variation of GSFE energy for a displacement u discretized on the hull
function gradient_Misfit(u::Array{Float64, 2}, hull::Hull)
    v = Reflection(u, hull) - u
    return (-hull.bl.iR) * gradient_GSFE(hull.G + hull.bl.R * v, hull.bl) - hull.bl.R * gradient_GSFE(hull.G + hull.bl.iR * v, hull.bl)
    # return 2 * hull.bl.R * gradient_GSFE(hull.G + v, hull.bl)
end

# Implements elastic energy for a displacement u discretized on the hull
function Elastic(u::Array{Float64, 2}, hull::Hull)
    N = hull.N
    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    @views return .5 / N^2 * real( dot(W[:], v[:]) )
end

# Implements first variation of elastic energy for a displacement u discretized on the hull
function gradient_Elastic(u::Array{Float64, 2}, hull::Hull)
    N = hull.N

    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    hull.iplan * W
    return reshape(real(W), (2, N^2))
end

# Computes the total (potential) energy of the bilayer for a displacement u discretized on the hull
function Energy(u::Array{Float64, 2}, hull::Hull)
    return 2*Elastic(u, hull) + Misfit(u, hull)
end

# Computes the first variation of the total (potential) energy of the bilayer for a displacement u discretized on the hull

function Gradient(u::Array{Float64, 2}, hull::Hull)
    grad = 2*gradient_Elastic(u, hull) + gradient_Misfit(u, hull)
    return grad
end

function Gradient!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    # N = hull.N
    #
    # v = hull.plan * reshape(u, (2, N, N))
    # W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    # hull.iplan * W
    # grad[:] = 2*real(W)
    #
    # # v = 2 * hull.bl.iR * u
    # # grad .+= 2 * hull.bl.R * gradient_GSFE(hull.G + v, hull.bl)
    #
    # v = Reflection(u, hull) - u
    # grad .+= (-hull.bl.iR) * gradient_GSFE(hull.G + hull.bl.R * v, hull.bl) - hull.bl.R * gradient_GSFE(hull.G + hull.bl.iR * v, hull.bl)
    grad .= gradient_Misfit(u, hull) + 2*gradient_Elastic(u, hull)
    nothing
end
