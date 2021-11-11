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
            -12.3222215429e-3 .* (sin.(s) .+ sin.(s+t)) .+
                -2.077e-3 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * -2) .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(s) .- cos.(s.+t)) .+
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*s)),
            -12.3222215429e-3 .* (sin.(t) .+ sin.(s+t)) .+
                -2.077e-3 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+
                (0.7834479467e-3 * -2) .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                2.3972862255e-3 .* (cos.(t) .- cos.(s.+t)) .+
                (0.25897164873e-3 * 2) .* (cos.(2*s .+ 2*t) .- cos.(2*t))
            );
end
