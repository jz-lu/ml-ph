### Implementation of GSFE for MoS2 0º homo-bilayers.
l = 3.122/sqrt(3)

# MoS2 elastic constants
K = 78.4478  # Bulk modulus
G = 63.0955  # Shear modulus

## NEW COEFFS (meV)
c0 = 3.702089895154718846e+01
c1 = 1.553573021770572815e+01
c2 = -4.358629720010219977e+00
c3 = -1.338304945182930306e+00
c4 = 0.000
c5 = 0.000

### Implementation of GSFE for MoS2 0º homo-bilayers.

l = 3.18 /sqrt(3)

# MoS2 elastic constants
K = 2.0*24.933  # Bulk modulus
G = 2.0*15.774  # Shear modulus

function GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]
    return c1 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .+
                c2 .* (cos.(s.+2. .*t) .+ cos.(s .-t) .+ cos.(2. .*s .+ t)) .+
                c3 .* (cos.(2. .*s) .+ cos.(2. .*t) .+ cos.(2. .*(s .+ t))) .+
                c4 .* (sin.(s) .+ sin.(t) .- sin.(s.+t)) .+ 
                c5 .* (sin.(2. .*s .+ 2. .*t) .- sin.(2. .*s) .- sin.(2. .*t))
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

    # return bl.invE' * vcat(
    #         -c1 .* (sin.(s) .+ sin.(s .+ t)) .+ 
    #             -c2 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+ 
    #             -2. .* c3 .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
    #             c4 .* (cos.(s) -. cos.(s .+ t)) .+
    #             2. .* c5 .* (cos.(2. .*(s .+ t)) .- cos.(2. .*s)), 
    #         -c1 .* (sin.(t) .+ sin.(s .+ t)) .+
    #             -c2 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+ 
    #             -2. .* c3 .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
    #             c4 .* (cos.(t) .- cos.(s .+ t)) .+
    #             2. .* c5 .* (cos.(2. .*(s .+ t)) .- cos.(2. .*t))
    #         );
end
