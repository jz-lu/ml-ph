# lattice constant, in angstrom
l = 1.42

# Graphene elastic constants (eV/unit cell)
K = 2.0*34.759  # Bulk modulus
G = 2.0*23.676  # Shear modulus

μ = G
λ = K - 2*μ/3

c1 = 4.063543396e-3
c2 = -0.374e-3
c3 = -0.094514083e-3
c4 = 0;
c5 = 0;

c_list = [c1, c2, c3, c4, c5]

function GSFE(X::Array{Float64,2}, bl::Bilayer)
    R = bl.invE*X
    s = R[1,:]
    t = R[2,:]

    return  c1 .* (cos.(s) .+ cos.(t) .+ cos.(s.+t)) .+
                c2 .* (cos.(s.+2. .*t) .+ cos.(s-t) .+ cos.(2. .*s .+ t)) .+
                c3 .* (cos.(2. .*s) .+ cos.(2. .*t) .+ cos.(2. .*(s .+ t)))
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
