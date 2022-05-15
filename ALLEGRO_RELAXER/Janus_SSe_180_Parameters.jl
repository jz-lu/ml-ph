## MoS2 + Janus TMDC 180 degree stack, S-Se interface

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


# GSFE coefficients
c1 = 13.39e-3
c2 = -2.323e-3
c3 = -0.703e-3
c4 = 3.726e-3
c5 = 0.4504e-3

c_list = [c1, c2, c3, c4, c5]

# height coefficients
d1 = 0.1185
d2 = -0.01608
d3 = -0.006753
d4 = 0.01927
d5 = 0.0004229

d_list = [d1, d2, d3, d4, d5]

κ = 0.21

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

# each G component of GSFE. c_m component, the n-th G component
function GSFE_comp(X::Array{Float64,2}, bl::Bilayer, m::Int64, n::Int64)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    if m <=3
        if m == 1
            n_list = [1 0;
                      0 1;
                      1 1];

        elseif m == 2
            n_list = [1 2;
                      1 -1;
                      2 1];

        elseif m == 3
            n_list = [2 0;
                      0 2;
                      2 2]
        end

        return c_list[m] .* cos.(n_list[n,1] .* s + n_list[n,2] .* t)
    else
        if m == 4
            n_list = [1 0;
                      0 1;
                      -1 -1];

        elseif m == 5
            n_list = [2 0;
                      0 2;
                      -2 -2]
        end
        return c_list[m] .* sin.(n_list[n,1] .* s + n_list[n,2] .* t)
    end
end



function inter_ph(bl::Bilayer, idx::Int64)
    Gvec = bl.iR * bl.invE'
    prefac = zeros(Float64, (2, 2, 3))

    if idx == 1 || idx == 4
        n_list = [1 0;
                  0 1;
                  -1 -1];

    elseif idx == 2
        n_list = [1 2;
                  1 -1;
                  2 1];

    elseif idx == 3 || idx == 5
        n_list = [2 0;
                  0 2;
                  -2 -2]
    end

    G = zeros(Float64, (2, 3))
    for i = 1:3
        G[:, i] = n_list[i, 1] * Gvec[:,1] + n_list[i, 2] * Gvec[:,2]
    end

    for a = 1:2, b = 1:2, c = 1:3
        prefac[a,b,c] = G[a,c] * G[b,c]
    end

    return prefac
end

function GSFE_height(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    return d1 .* (cos.(s) .+ cos.(t) .+ cos.(s .+ t)) .+
                d2 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .+
                d3 .* (cos.(2*s) .+ cos.(2*t) .+ cos.(2. .*(s .+ t))) .+
                d4 .* (sin.(s) .+ sin.(t) .- sin.(s.+t)) .+
                d5 .* (sin.(2*s .+ 2*t) .- sin.(2*s) .- sin.(2*t))
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

function gradient_GSFE_height(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]'
    t = R[2,:]'

    return bl.invE' * vcat(
            -d1 .* (sin.(s) .+ sin.(s+t)) .+
                 -d2 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
                (-d3 * 2) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                d4 .* (cos.(s) .- cos.(s.+t)) .+
                (d5 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
            -d1 .* (sin.(t) .+ sin.(s+t)) .+
                -d2 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
                (-d3 * 2) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                d4 .* (cos.(t) .- cos.(s.+t)) .+
                (d5 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
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
