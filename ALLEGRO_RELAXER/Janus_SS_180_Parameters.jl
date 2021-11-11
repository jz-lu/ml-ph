## MoS2 + Janus TMDC, S-S interface

# lattice constant, in angstrom
l = 3.1 /sqrt(3)

# MoSSe
G1 = 28.5367
K1 = 49.035

# MoS2
G2 = 32.8647
K2 = 50.4193

K = 0.5 * (K1 + K2)
G = 0.5 * (G1 + G2)

λ = K - 2*G/3
μ = G

c1 = 15.68e-3
c2 = -3.012e-3
c3 = -0.8563e-3
c4 = 0.996e-3
c5 = 0.1901e-3

c_list = [c1, c2, c3, c4, c5]

function GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    return c1 .* (cos.(s) .+ cos.(t) .+ cos.(s .+ t)) .+
                c2 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .+
                c3 .* (cos.(2*s) .+ cos.(2*t) .+ cos.(2. .*(s .+ t))) .+
                c4 .* (sin.(s) .+ sin.(t) .- sin.(s.+t)) .+
                c5 .* (sin.(2*s .+ 2*t) .- sin.(2*s) .- sin.(2*t))
end


function gradient_GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]'
    t = R[2,:]'

    return bl.invE' * vcat(
            -c1 .* (sin.(s) .+ sin.(s+t)) .+
                 -c2 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
                (-c3 * 2) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                c4 .* (cos.(s) .- cos.(s.+t)) .+
                (c5 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
            -c1 .* (sin.(t) .+ sin.(s+t)) .+
                -c2 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
                (-c3 * 2) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                c4 .* (cos.(t) .- cos.(s.+t)) .+
                (c5 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
            );
end

function Gradient!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    N = hull.N

    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    hull.iplan * W
    grad[:] = 2*real(W)

    v = 2 * hull.bl.iR * u

    grad .+= 2 .* hull.bl.R * gradient_GSFE(hull.G + v, hull.bl)
    nothing
end
