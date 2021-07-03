### Implementation of GSFE for MoS2 0º homo-bilayers.

# MoS2 elastic constants
K = 2.0*24.933  # Bulk modulus
G = 2.0*15.774  # Shear modulus

function GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    return 14.020e-3 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .-
                2.542e-3 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .- 
                0.884e-3 .* (cos.(2*s) .+ cos.(2*t) .+ cos.(2. .*(s .+ t)))
end

# MoS2 GSFE functional gradient
function gradient_GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]'
    t = R[2,:]'

    return bl.invE' * vcat(
            -14.020e-3 .* (sin.(s) .+ sin.(s+t)) .+ 
                2.542e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+ 
                (0.884e-3 * 2.) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))),
            -14.020e-3 .* (sin.(t) .+ sin.(s+t)) .+
                2.542e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+ 
                (0.884e-3 * 2.) .* (sin.(2*t) .+ sin.(2. .*(s .+ t)))
            );
end

## Efficient in-place implementation of the total energy gradient.

function Gradient!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    N = hull.N

    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    hull.iplan * W
    grad[:] = 2*real(W) 

    S = Reflection(u, hull) - u
    R = hull.Gref + (hull.bl.invE*hull.bl.R) * S
    @views s = R[1,:]'
    @views t = R[2,:]'
    grad .+= (hull.bl.iR * hull.bl.invE') * vcat(
            14.020e-3 .* (sin.(s) .+ sin.(s.+t)) .-
                2.542e-3 .* (sin.(s.+(2.).*t) .+ sin.(s.-t) .+ 2. .*sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*s) .+ sin.((2.).*(s.+t))),
            14.020e-3 .* (sin.(t) .+ sin.(s.+t)) .-
                2.542e-3 .* (2. .*sin.(s.+(2.)*t) .+ sin.(t.-s) .+ sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*t) .+ sin.((2.).*(s.+t)))
            )

    R = hull.Gref + (hull.bl.invE*hull.bl.iR) * S
    @views s = R[1,:]'
    @views t = R[2,:]'
    grad .+= (hull.bl.R * hull.bl.invE') * vcat(
            14.020e-3 .* (sin.(s) .+ sin.(s.+t)) .-
                2.542e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*s) .+ sin.((2.).*(s.+t))),
            14.020e-3 .* (sin.(t) .+ sin.(s.+t)) .-
                2.542e-3 .* ((2.).*sin.(s.+(2.).*t) .+ sin.(t.-s) .+ sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*t) .+ sin.((2.).*(s.+t)))
            )
    nothing
end

function Gradient_nosym!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    N = hull.N

    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    hull.iplan * W
    grad[:] = 2*real(W) 

    R = hull.Gref - (2 * hull.bl.invE*hull.bl.iR) * u
    @views s = R[1,:]'
    @views t = R[2,:]'
    grad .+= (2 * hull.bl.R * hull.bl.invE') * vcat(
            14.020e-3 .* (sin.(s) .+ sin.(s.+t)) .-
                2.542e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*s) .+ sin.((2.).*(s.+t))),
            14.020e-3 .* (sin.(t) .+ sin.(s.+t)) .-
                2.542e-3 .* ((2.).*sin.(s.+(2.).*t) .+ sin.(t.-s) .+ sin.((2.).*s.+t)) .-
                (0.884e-3 * 2.) .* (sin.((2.).*t) .+ sin.((2.).*(s.+t)))
            )
    nothing
end

