### Implementation of GSFE for MoS2 0º homo-bilayers.

# lattice constant, in angstrom
l = 3.18 /sqrt(3)

# MoS2 elastic constants
K = 2.0*24.933  # Bulk modulus
G = 2.0*15.774  # Shear modulus


# MoS2 GSFE functional
function GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    return 12.3222215429e-3 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .+
                -2.077e-3 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .+
                -0.7834479467e-3 .* (cos.(2*s) .+ cos.(2*t) .+ cos.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (sin.(s) .+ sin.(t) .- sin.(s.+t)) .+
                0.25897164873e-3 .* (sin.(2*s .+ 2*t) .- sin.(2*s) .- sin.(2*t))
end

# MoS2 GSFE functional gradient
function gradient_GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]'
    t = R[2,:]'

    return bl.invE' * vcat(
        # ! partial F/partial s
            -12.3222215429e-3 .* (sin.(s) .+ sin.(s+t)) .+
                 2.077e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * 2) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(s) .- cos.(s.+t)) .+ # TODO this and next line
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
            # ! partial F/partial t
            -12.3222215429e-3 .* (sin.(t) .+ sin.(s+t)) .+
                2.077e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * 2) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(t) .- cos.(s.+t)) .+
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
            );
end

## Efficient in-place implementation of the total energy gradient.

# function Gradient!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
#     N = hull.N
#
#     v = hull.plan * reshape(u, (2, N, N))
#     W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
#     hull.iplan * W
#     grad[:] = 2*real(W)
#
#     S = Reflection(u, hull) - u
#     R = hull.Gref + (hull.bl.invE*hull.bl.R) * S
#     @views s = R[1,:]'
#     @views t = R[2,:]'
#     grad .+= (-hull.bl.iR * hull.bl.invE') * vcat(
#             -12.3222215429e-3 .* (sin.(s) .+ sin.(s+t)) .+
#                  2.077e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
#                 (0.7834479467e-3 * 2.) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
#                 2.3972862255e-3 .* (cos.(s) .- cos.(s.+t)) .+
#                 (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
#             -12.3222215429e-3 .* (sin.(t) .+ sin.(s+t)) .+
#                 2.077e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
#                 (0.7834479467e-3 * 2.) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
#                 2.3972862255e-3 .* (cos.(t) .- cos.(s.+t)) .+
#                 (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
#             )
#
#     R = hull.Gref + (hull.bl.invE*hull.bl.iR) * S
#     @views s = R[1,:]'
#     @views t = R[2,:]'
#     grad .+= (-hull.bl.R * hull.bl.invE') * vcat(
#             -12.3222215429e-3 .* (sin.(s) .+ sin.(s+t)) .+
#                  2.077e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
#                 (0.7834479467e-3 * 2.) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
#                 2.3972862255e-3 .* (cos.(s) .- cos.(s.+t)) .+
#                 (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
#             -12.3222215429e-3 .* (sin.(t) .+ sin.(s+t)) .+
#                 2.077e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
#                 (0.7834479467e-3 * 2.) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
#                 2.3972862255e-3 .* (cos.(t) .- cos.(s.+t)) .+
#                 (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
#             )
#     nothing
# end

function Gradient_nosym!(grad::Array{Float64, 2}, u::Array{Float64, 2}, hull::Hull)
    N = hull.N

    v = hull.plan * reshape(u, (2, N, N))
    W = dropdims(sum(hull.Hess_Elastic .* reshape(v, (2, 1, N, N)), dims=1), dims=1);
    hull.iplan * W
    grad[:] = 2*real(W)

    R = hull.Gref - (2 * hull.bl.invE*hull.bl.iR) * u
    @views s = R[1,:]'
    @views t = R[2,:]'
    grad .+= (-2 * hull.bl.R * hull.bl.invE') * vcat(
            -12.3222215429e-3 .* (sin.(s) .+ sin.(s+t)) .+
                 2.077e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * 2.) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(s) .- cos.(s.+t)) .+
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
            -12.3222215429e-3 .* (sin.(t) .+ sin.(s+t)) .+
                2.077e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * 2.) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(t) .- cos.(s.+t)) .+
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
            )
    nothing
end
