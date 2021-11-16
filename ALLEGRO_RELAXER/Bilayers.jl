### Relaxation of bilayers using Hull framework.
## Common setup to all rotated homo-bilayers.

using LinearAlgebra

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
    tC1::Array{Float64,4}
    tC2::Array{Float64,4}

    function Bilayer(l, θ, K1, G1, K2, G2)
        # E = l * [3/2 3/2; -sqrt(3)/2 sqrt(3)/2]
        E = l * [1 1/2 ; 0 sqrt(3)/2] * sqrt(3)
        P = l * [1; 0]

        ## Symmetries and rotations:
        R  = [cos(θ/2) sin(θ/2); -sin(θ/2) cos(θ/2)]
        iR = [cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)]

        ## Lattice parameters for twisted layers:
        tE = [iR*E,R*E]
        tP = [iR*P,R*P]

        ## Intra-layer elastic tensor, with stress tensor linear ordering [E_11; E_12; E_21; E_22]
        C1 = [K1+G1 0 0 K1-G1; 0 G1 G1 0; 0 G1 G1 0; K1-G1 0 0 K1+G1]
        C2 = [K2+G2 0 0 K2-G2; 0 G2 G2 0; 0 G2 G2 0; K2-G2 0 0 K2+K2]
        T = inv(E)  * (R - iR)
        T = [T[1,1] 0 T[1,2] 0; 0 T[1,1] 0 T[1,2]; T[2,1] 0 T[2,2] 0; 0 T[2,1] 0 T[2,2]] # Hull [0,1]^2 -> Real space gradient transformation
        tC1 =  T * C1 * T';
        tC2 =  T * C2 * T';
        tC1 = reshape(tC1, (2,2,2,2))
        tC2 = reshape(tC2, (2,2,2,2))

        new(l, E, 2*pi*inv(E), P, θ, R, iR, tE, tP, tC1, tC2)
    end
end

# Discretized Hull structure
mutable struct Hull
    bl::Bilayer
    N::Int
    G::Array{Float64,2}
    Gref::Array{Float64,2}
    Hess_Elastic1::Array{Float64,4}
    Hess_Elastic2::Array{Float64,4}
    plan
    iplan

    # Generate Hull discretization
    # modified Hull: offset Hess_Elastic by the center site index in the k-space direct coordinate
    function Hull(bl::Bilayer, N, kidx::Array{Float64,1})
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

        Hess1 = zeros(2,2,N,N) # layer 1
        Hess2 = zeros(2,2,N,N) # layer 2
        K1 = K1 .+ 2*pi*kidx[1]
        K2 = K2 .+ 2*pi*kidx[2]
        for a=1:2, c=1:2
            Hess1[a,c,:,:] =  bl.tC1[a,1,c,1] * (K1.*K1) .+ (bl.tC1[a,2,c,1] .+ 
                bl.tC1[a,1,c,2]) * (K1.*K2) .+ bl.tC1[a,2,c,2] * (K2.*K2)
            Hess2[a,c,:,:] =  bl.tC2[a,1,c,1] * (K1.*K1) .+ (bl.tC2[a,2,c,1] .+ 
                bl.tC2[a,1,c,2]) * (K1.*K2) .+ bl.tC2[a,2,c,2] * (K2.*K2)
        end

        new(bl, N, G, Gref, Hess1, Hess2, plan, iplan)
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
    u1 = u[1:2, :]
    u2 = u[3:4, :]
    return sum(GSFE(hull.G + hull.bl.iR * (u1+u2), hull.bl))
end


# Implements first variation of GSFE energy for a displacement u discretized on the hull
function gradient_Misfit(u::Array{Float64, 2}, hull::Hull)
    u1 = u[1:2, :]
    u2 = u[3:4, :]
    v = hull.bl.iR * (u1+u2)
    grad1 = hull.bl.R * gradient_GSFE(hull.G + v, hull.bl)
    return [grad1; grad1]
end

# Implements elastic energy for a displacement u discretized on the hull
function Elastic(u::Array{Float64, 2}, hull::Hull)
    u1 = u[1:2, :]
    u2 = u[3:4, :]
    N = hull.N
    # v = hull.plan * reshape(u, (2, N, N))
    v1 = hull.plan * reshape(u1, (2, N, N))
    v2 = hull.plan * reshape(u2, (2, N, N))
    # W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    W1 = dropdims(sum(hull.Hess_Elastic1 .* reshape(v1, (2, 1, N, N)), dims=1), dims=1);
    W2 = dropdims(sum(hull.Hess_Elastic2 .* reshape(v2, (2, 1, N, N)), dims=1), dims=1);
    # @views return .5 / N^2 * real( dot(W[:], v[:]) )
    @views return .5 / N^2 * (real( dot(W1[:], v1[:]) ) + real( dot(W2[:], v2[:]) ))
end

# Implements first variation of elastic energy for a displacement u discretized on the hull
function gradient_Elastic(u::Array{Float64, 2}, hull::Hull)
    u1 = u[1:2, :]
    u2 = u[3:4, :]
    N = hull.N
    # v = hull.plan * reshape(u, (2, N, N))
    v1 = hull.plan * reshape(u1, (2, N, N))
    v2 = hull.plan * reshape(u2, (2, N, N))
    # W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    W1 = dropdims(sum(hull.Hess_Elastic1 .* reshape(v1, (2, 1, N, N)), dims=1), dims=1);
    W2 = dropdims(sum(hull.Hess_Elastic2 .* reshape(v2, (2, 1, N, N)), dims=1), dims=1);
    # hull.iplan * W
    hull.iplan * W1
    hull.iplan * W2
    # return reshape(real(W), (2, N^2))
    return reshape(real([W1; W2]), (4, N^2))
end

# Computes the total (potential) energy of the bilayer for a displacement u discretized on the hull
function Energy(u::Array{Float64, 2}, hull::Hull)
    return Elastic(u, hull) + Misfit(u, hull)
end

# Computes the first variation of the total (potential) energy of the bilayer for a displacement u discretized on the hull
function Gradient(u::Array{Float64, 2}, hull::Hull)
    return gradient_Elastic(u, hull) + gradient_Misfit(u, hull)
end


function Gradient!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    grad .= gradient_Misfit(u, hull) + gradient_Elastic(u, hull)
    nothing
end
