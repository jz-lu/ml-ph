using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
include("gsfe_func.jl")

using PyPlot
using PyCall
using Interpolations
using Optim
using LinearAlgebra
using ArgParse
using NPZ
using Base.Iterators

# Obtain the DOS Gaussian standard deviation based on twist angle
function DOS_stdev(θdeg)
    width = 0.45
    if θdeg < 1.01
        width = 0.03
    elseif θdeg < 1.51
        width = 0.09
    elseif θdeg < 2.01
        width = 0.13
    elseif θdeg < 2.51
        width = 0.15
    elseif θdeg < 3.51
        width = 0.21
    elseif θdeg < 5.51
        width = 0.27
    elseif θdeg < 6.51
        width = 0.31
    elseif θdeg < 7.51
        width = 0.37
    elseif θdeg < 9.01
        width = 0.43
    end
    return width / (2 * sqrt(2 * log(2))) / 5
end

partition_density = 1000
function compute_DOS(weights, θdeg, omegas_k, A=1)
    println("Computing DOS... ")
    variance = DOS_stdev(θdeg)^2
    omegas = LinRange(0, 3, partition_density)
    return omegas, [sum(
        sum(
            [weight * A * exp(-(omega - omega_k)^2 / (2 * variance)) 
            for (omega_k, weight) in zip(omegas_k, weights) ]
        )
    ) for omega in omegas]
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--deg", "-d"
            help = "twist angle in degrees"
            # required = true
            arg_type = Float64
            default = 2.
        "-N", "-n"
            help = "grid size (N if grid is NxN)"
            arg_type = Int
            default = 51
        "--mesh", "-m"
            help = "DOS over mesh instead of bands over line"
            action = :store_true
        "--cut", "-c"
            help = "Grid shell cutoff"
            arg_type = Int
            default = 6
        "--out", "-o"
            help = "output directory path (must end in slash)"
            arg_type = String
            default = "."
    end

    return parse_args(s)
end

close("all")

## Parameters:
parsed_args = parse_commandline()
dir = string(abspath(parsed_args["out"]), "/")
θ = deg2rad(parsed_args["deg"])
θdeg = parsed_args["deg"]
N = parsed_args["N"]
hN = div(N,2)
blg = Bilayer(l, θ, K, G)
hull = Hull(blg, N, [0.,0.]);

# Cutofff radius in k-space
G_cutoff = parsed_args["cut"] # need to be large for small angles
# number of points to sample on the high symmetry line
nq = 100
nr = 50
mode = "line"
if parsed_args["mesh"] == true
    mode = "mesh"
    nq = 45
end
wf_idx = 1

## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
f(u) = Energy(u, hull)
g!(storage, u) = Gradient!(storage, u, hull)
@time results = optimize(f, g!, zeros(2, N^2), LBFGS(), Optim.Options(allow_f_increases=true))

# static solution
u0 = reshape(Optim.minimizer(results), (2, N^2))
print(results)
# u0 =  blg.E * [1/3; 1/3 ]
# ∇u = ∇(u0, hull)
# Fourier component of u
uG = hull.plan * reshape(u0, (2,N,N))


# moire cell
EM = inv(blg.R-blg.iR) * blg.E
# moire reciprocal cell
GM = 2*pi*inv(EM)'
GM1 = GM[:,1]
GM2 = GM[:,2]
# generate frequency grid
freq = fftfreq(N,N) # the second argument: sampling rate of the input argument = a0*N
ind = sortperm(freq)
uG = uG[:,ind,:]
uG = uG[:,:,ind]
freq = sort(freq)
f1 = freq .* ones(1,N)
f2 = freq' .* ones(N,1)
k =  GM * [f1[:]'; f2[:]']
kx = reshape(k[1,:], (N,N))
ky = reshape(k[2,:], (N,N))

# define mesh grid
xgrid = reshape(range(0,stop=nq-1,step=1), (1,nq))
ygrid = reshape(range(0,stop=nq-1,step=1), (nq,1))
xmesh = xgrid .* ones(nq,1)
ymesh = ygrid .* ones(1,nq)

# high symmetry points
if mode == "line"
    Gpt = [0, 0]
    Kpt = (GM1[1:2]+2*GM2[1:2])/3
    Mpt = 1/2*(GM1[1:2] + GM2[1:2])

    Kpt_idx = [1, 2]/3
    Gpt_idx = [0, 0]
    Mpt_idx = [1/2, 1/2] # indices of monolayer high symmetry point (not used)

    l1 = zeros(2,nq+1)
    l2 = zeros(2,nq+1)
    l3 = zeros(2,nq)
    i1 = zeros(2,nq+1)
    i2 = zeros(2,nq+1)
    i3 = zeros(2,nq)
    kline = zeros(2,nq*3)
    kline_idx = zeros(2,nq*3)
    for j = 1:2
        l1[j,:] = range(Kpt[j],Gpt[j],length=nq+1)
        l2[j,:] = range(Gpt[j],Mpt[j],length=nq+1)
        l3[j,:] = range(Mpt[j],Kpt[j],length=nq)
        i1[j,:] = range(Kpt_idx[j],Gpt_idx[j],length=nq+1)
        i2[j,:] = range(Gpt_idx[j],Mpt_idx[j],length=nq+1)
        i3[j,:] = range(Mpt_idx[j],Kpt_idx[j],length=nq)
        kline[j,:] = [l1[j,1:nq];l2[j,1:nq];l3[j,:]]
        kline_idx[j,:] = [i1[j,1:nq];i2[j,1:nq];i3[j,:]]
    end

elseif mode == "mesh"
    # create a grid in k space
    kline = GM * [xmesh[:]' .- max(xmesh[:]...)/2;ymesh[:]' .- max(ymesh[:]...)/2]/nq
    kline_idx = [xmesh[:]' .- max(xmesh[:]...)/2;ymesh[:]' .- max(ymesh[:]...)/2]/nq

    kline = GM * [xmesh[:]';ymesh[:]']/nq
    kline_idx = [xmesh[:]';ymesh[:]']/nq
end

# create a grid of reciprocal lattice vectors within the cutoff radius
Gcut = G_cutoff * norm(GM[:,1])
knorm = sqrt.(k[1,:].^2+k[2,:].^2)
id = findall(x -> x < Gcut*1.1, knorm)
q_list = k[:,id]

KG = reshape(hull.Hess_Elastic,(2,2,N^2))
Kq = KG[:,:,id]

# find the indices of the center site
id_q1 = findall(x -> x == 0, q_list[1,:])
id_q2 = findall(x -> x == 0, q_list[2,:])
tmp_id = findall(x -> x in id_q2, id_q1)
q0_idx = id_q1[tmp_id][1]

##
evals = zeros(Float64,(size(kline,2),size(q_list,2)*2))
# interlayer
# first calculate the Fourier component of the gradient of the relaxed misfit energy
A = zeros(Float64, (N, N))
plan = plan_bfft(A, (1,2) )

# Fourier decomposition of relaxed GSFE
B = hull.G + 2 * blg.iR * u0 # relaxed configuration space position
VG = zeros(ComplexF64, (2, 2, N, N))

for ci = 1:5
    prefac_tmp = inter_ph(blg, ci)
    gsfe = zeros(ComplexF64, (N, N, 3))

    for j = 1:3
        gsfe[:,:,j] = plan * reshape(GSFE_comp(B, blg, ci, j), (N, N))
    end

    gsfe=gsfe[ind,:,:]
    gsfe=gsfe[:,ind,:] # sort Fourier component of GSFE by G vector
    # pcolormesh(kx,ky,real.(gsfe[:,:,1]+gsfe[:,:,2]+gsfe[:,:,3]))

    for c = 1:3
        for a = 1:2, b = 1:2
            global VG[a,b,:,:] = VG[a,b,:,:] - 2 * prefac_tmp[a,b,c] * gsfe[:,:,c] / N^2
        end
    end
end

VG = reshape(VG, (2,2,N^2))

Hinter = zeros(ComplexF64,(size(q_list,2)*2, size(q_list,2)*2))
qidx = zeros(Int64, (size(q_list,2), size(q_list,2)))
dG = zeros(Float64, (2, size(q_list,2), size(q_list,2)))

tol = 1e-6
for q1_idx = 1:size(q_list,2)
    for q2_idx = 1:size(q_list,2)

        id1 = q1_idx*2-1;
        id2 = q2_idx*2-1;
        dq = q_list[:,q2_idx]-q_list[:,q1_idx]
        dG[:, q1_idx, q2_idx] = dq
        id_q1 = findall(x -> abs(x-dq[1])<tol, k[1,:])
        id_q2 = findall(x -> abs(x-dq[2])<tol, k[2,:])

        tmp_id = findall(x -> x in id_q2, id_q1)
        if !isempty(tmp_id)
            qid = id_q1[tmp_id]
            qidx[q1_idx, q2_idx] = qid[1]

            Hinter[id1:id1+1,id2:id2+1]=VG[:,:,qid]

        end
    end
end

a0 = l*sqrt(3)
A0 = sqrt(3)/2*a0^2

evecs_all = zeros(ComplexF64,(size(q_list,2)*2, size(q_list,2)*2, size(kline,2)))
ievecs_all = zeros(ComplexF64,(size(q_list,2)*2, size(q_list,2)*2, size(kline,2)))
ievals = zeros(Float64,(size(kline,2),size(q_list,2)*2))

for k_idx = 1:size(kline,2)
    Hintra = zeros(ComplexF64,size(q_list,2)*2, size(q_list,2)*2)

    KG = zeros(Float64, (2,2,size(q_list,2)))
    for g_idx = 1:size(KG, 3)
        Gx = q_list[1,g_idx] + kline[1, k_idx]
        Gy = q_list[2,g_idx] + kline[2, k_idx]

        KG[:,:,g_idx] = [(λ+2*μ)*Gx^2+μ*Gy^2  (λ+μ)*Gx*Gy
                        (λ+μ)*Gx*Gy (λ+2*μ)*Gy^2+μ*Gx^2]

        id1 = g_idx*2-1
        id2 = id1+1
        Hintra[id1:id2,id1:id2]=KG[:,:,g_idx]

    end

    global H = (Hinter+Hintra)/A0 # convert from eV / cell to eV / A^2

    # intralayer part
    iF = eigen(Hintra)
    ieid = sortperm(real.(iF.values))
    ievals[k_idx,:] = real(iF.values[ieid])
    ievecs_all[:, :, k_idx] = iF.vectors[:,ieid]

    # the whole thing
    global F = eigen(H) # default sort by real(λ)
    eid = sortperm(real.(F.values))
    evals[k_idx,:] = real(F.values[eid])
    global evecs_all[:, :, k_idx] = F.vectors[:,eid]
    # The kth eigenvector can be obtained from the slice F.vectors[:, k].

    if mod(k_idx,10) == 0
        println("Diagonalization done ", k_idx, "/", size(kline,2))
    end
end

# normalization
rho = 7.61 * 1e-7 # area density of single-layer graphene in kg/m^2

w = sqrt.(abs.(evals) * 1e20 * 16.02   ./ rho) # in Hz
iw = sqrt.(abs.(ievals) * 1e20 * 16.02   ./ rho) # in Hz
lambda = λ/A0# in eV / A^2
lambda = lambda * 1.602e-19 *1e20; # in J / m^2
LM = a0 / θ *1e-10 # moire length in m
w0 = sqrt(lambda/rho) * (2*pi) / LM # in 1/s
# convert to meV
hbar = 4.1357e-15 # in eV.s
w0_meV = w0 * hbar *1e3
w0_THz = w0 / 10e12

if mode == "line"
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 15

    fig=figure(2, figsize = (6, 8))
    clf()
    xlabel("k")
    ylabel("ω/ω₀")
    plot(w/w0,color="black")
    title(string("θ = ", θdeg))
    ylim((0.0,3))
    xlim((0,size(w)[1]))
    figure(2, figsize = (6,8))
    savefig(string(dir, "band_koshino_", θdeg, ".png"))
    npzwrite(string(dir, "band_koshino_", θdeg, ".npz"), w/w0)
    println(string("Saved to ", string(dir, "band_koshino_", θdeg, ".png")))
elseif mode == "mesh"
    ifreqs = iw/w0
    iG0_piece = ievecs_all[2*q0_idx-1 : 2*q0_idx, :, :]
    iweights = transpose(dropdims(sum(abs.(iG0_piece).^2, dims=1), dims=1))
    println("iWeights size:")
    println(size(iweights))
    println("iFreqs size:")
    println(size(ifreqs))
    println("Flattening...")
    ifreqs = ifreqs[:]
    iweights = iweights[:]
    println("iWeights size (flat):")
    println(size(iweights))
    println("iFreqs size (flat):")
    println(size(ifreqs))
    iomegas, iDOS = compute_DOS(iweights, θdeg, ifreqs, 1)
    normalizer = maximum(iDOS)
    println("Intralayer done.")

    freqs = w/w0
    G0_piece = evecs_all[2*q0_idx-1 : 2*q0_idx, :, :]
    weights = transpose(dropdims(sum(abs.(G0_piece).^2, dims=1), dims=1))
    println("Weights size:")
    println(size(weights))
    println("iFreqs size:")
    println(size(freqs))
    freqs = freqs[:]
    weights = weights[:]
    println("Weights size (flat):")
    println(size(weights))
    println("Freqs size (flat):")
    println(size(freqs))
    omegas, DOS = compute_DOS(weights, θdeg, freqs, normalizer)

    # Plot the DOS
    println("Plotting DOS...")
    clf()
    figure(2, figsize = (6,8))
    xlabel("DOS")
    ylim((0.0, 3))
    plot(DOS, omegas, color="black")
    title(string("θ = ", θdeg))
    savefig(string(dir, "dos_koshino_", θdeg, ".png"))
    npzwrite(string(dir, "dos_koshino_", θdeg, ".npz"), DOS)
    npzwrite(string(dir, "omegas_koshino_", θdeg, ".npz"), omegas)
    println(string("Saved to ", string(dir, "dos_koshino_", θdeg, ".png")))
end