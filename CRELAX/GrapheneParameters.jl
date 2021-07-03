### Implementation of GSFE for graphene homo-bilayers.

# lattice constant, in angstrom
l = 1.42

# Graphene elastic constants (eV/unit cell)
K = 2*34.759  # Bulk modulus
G = 2*23.676  # Shear modulus

# From Koshino 2017
# K1 = 40.8626
# G1 = 45.5776
#
# K2 = 40.8626
# G2 = 45.5776

c1 = 4.063543396e-3
c2 = -0.374e-3
c3 = -0.094514083e-3
c4 = 0;
c5 = 0;

# Graphene GSFE functional
# input: real space positions
function GSFE(X::Array{Float64,2}, bl::Bilayer)

    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]

    return  c1 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .+
                c2 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .+
                c3 .* (cos.(2. .*s) .+ cos.(2. .*t) .+ cos.(2. .*(s .+ t)))
end

function gradient_GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]

    G = similar(X)
    for i=1:size(X,2)
      @inbounds G[1,i] = -c1 * (sin(s[i]) + sin(s[i]+t[i])) +
                 -c2 * (sin(s[i]+2. *t[i]) + sin(s[i]-t[i]) + 2. *sin(2. *s[i] + t[i])) +
                 -c3 * (sin(2. *s[i]) + sin(2. *(s[i]+t[i])))
      @inbounds G[2,i] = -c1 * (sin(t[i]) + sin(s[i]+t[i])) +
                  -c2 * (2. *sin(s[i]+2. *t[i]) + sin(t[i]-s[i]) + sin(2. *s[i]+t[i])) +
                  -c3 * (sin(2. *t[i]) + sin(2. *(s[i]+t[i])))
    end

    return (bl.invE') * G
end

#
