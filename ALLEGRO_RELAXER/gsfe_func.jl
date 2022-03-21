include("MoSe2_WSe2_0_Parameters.jl")

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
            -c1 .* (sin.(s) .+ sin.(s .+ t)) .+ 
                -c2 .* (sin.(s.+2. .*t) .+ sin.(s.-t) .+ 2. .*sin.(2. .*s .+ t)) .+ 
                -2. .* c3 .* (sin.(2*s) .+ sin.(2. .*(s .+ t))) .+
                c4 .* (cos.(s) .- cos.(s .+ t)) .+
                2. .* c5 .* (cos.(2. .*(s .+ t)) .- cos.(2. .*s)), 
            -c1 .* (sin.(t) .+ sin.(s .+ t)) .+
                -c2 .* (2*sin.(s.+2. .*t) .+ sin.(t.-s) .+ sin.(2. .*s .+ t)) .+ 
                -2. .* c3 .* (sin.(2*t) .+ sin.(2. .*(s .+ t))) .+
                c4 .* (cos.(t) .- cos.(s .+ t)) .+
                2. .* c5 .* (cos.(2. .*(s .+ t)) .- cos.(2. .*t))
            );
end
