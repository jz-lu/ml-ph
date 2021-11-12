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

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--deg", "-d"
            help = "twist angle in degrees"
            required = true
            arg_type = Float64
        "-N", "-n"
            help = "grid size (N if grid is NxN)"
            arg_type = Int
            required = true
        "--cut", "-c"
            help = "Grid shell cutoff"
            arg_type = Int
            default = 7
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
θ = deg2rad(parsed_args["deg"])
N = parsed_args["N"]
hN = div(N,2)
blg = Bilayer(l, θ, K, G)
hull = Hull(blg, N, [0.,0.]);

# Cutofff radius in k-space
G_cutoff = parsed_args["cut"] # need to be large for small angles
# number of points to sample on the high symmetry line
nq = 20
nr = 50
mode = "line" # mode = "mesh" or "line", "line" to calculate band structures, "mesh" for DOS (not implemented)
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

fig = figure(1,figsize = (8, 10))
fig.subplots_adjust(hspace=0.1, wspace=0.1)
clf()
subplot(2,1,1)
pcolormesh(kx,ky,abs.(uG[1,:,:]))
xlim((-5,5))
ylim((-5,5))
ax=gca()
ax.set_aspect(1)
colorbar()

subplot(2,1,2)
pcolormesh(kx,ky,abs.(uG[2,:,:]))
xlim((-5,5))
ylim((-5,5))
ax=gca()
ax.set_aspect(1)
colorbar()

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

    figure(10)
    # scatter(Kpt[1], Kpt[2])
    # scatter(Mpt[1], Mpt[2])
    # scatter(Gpt[1], Gpt[2])
    # scatter(K1[1], K1[2])
    # scatter(K2[1], K2[2])
    # scatter(G0[1,1], G0[2,1])
    # scatter(G0[1,2], G0[2,2])
    scatter(kline[1,:], kline[2,:])
    ax=gca()
    ax.set_aspect(1)
    figure(10)


elseif mode == "mesh"
    # create a grid in k space

    # kline = G0 * [xmesh[:]'; ymesh[:]'] / nq
    # kline_idx = [xmesh[:]'; ymesh[:]'] / nq
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

# no cutoff
# q_list = k
# Kq = KG

##
evals = zeros(Float64,(size(kline,2),size(q_list,2)*2))
# interlayer
# first calculate the Fourier component of the gradient of the relaxed misfit energy
A = zeros(Float64, (N, N))
plan = plan_bfft(A, (1,2) )

# Fourier decomposition of relaxed GSFE

B = hull.G + 2 * blg.iR * u0 # relaxed configuration space position
VG = zeros(ComplexF64, (2, 2, N, N))

# figure(3)

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

# ax=gca()
# ax.set_aspect(1)
# colorbar()
# figure(3)


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
            # println(qid)
            # println([q1_idx, q2_idx], dq, qid)
            # println(qid)
            qidx[q1_idx, q2_idx] = qid[1]

            Hinter[id1:id1+1,id2:id2+1]=VG[:,:,qid]

        end
    end
end


a0 = l*sqrt(3)
A0 = sqrt(3)/2*a0^2

evecs_all = zeros(ComplexF64,(size(q_list,2)*2, size(q_list,2)*2, size(kline,2)))

for k_idx = 1:size(kline,2)
    Hintra = zeros(ComplexF64,size(q_list,2)*2, size(q_list,2)*2)

    # an equivalent way to create KG
    # hull = Hull(blg, N, kline_idx[:,k_idx])
    # KG = hull.Hess_Elastic
    # KG = reshape(KG, (2, 2, N, N))
    # KG = KG[:,:,ind,:]
    # KG = KG[:,:,:,ind] # sort KG from smallest G vector to the largest
    # KG = reshape(KG,(2,2,N^2))
    # KG = KG[:,:,id] # within the cutoff radius

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

    # global H = Hintra/A0

    global F = eigen(H) # default sort by (real(λ), real(λ))
    eid = sortperm(real.(F.values))
    evals[k_idx,:] = real(F.values[eid])
    global evecs_all[:, :, k_idx] = F.vectors[:,eid]
    # The kth eigenvector can be obtained from the slice F.vectors[:, k].

    if mod(k_idx,10) == 0
        println("Diagonalizatsion done ", k_idx, "/", size(kline,2))
    end
end

# normalization
rho = 7.61 * 1e-7 # area density of single-layer graphene in kg/m^2
evals_shift = evals .- min(evals[:]...)
evals_shift = evals

w = sqrt.(abs.(evals_shift) * 1e20 * 16.02   ./ rho) # in Hz
lambda = λ/A0# in eV / A^2
lambda = lambda * 1.602e-19 *1e20; # in J / m^2
LM = a0 / θ *1e-10 # moire length in m
w0 = sqrt(lambda/rho) * (2*pi) / LM # in 1/s
# convert to meV
hbar = 4.1357e-15 # in eV.s
w0_meV = w0 * hbar *1e3
w0_THz = w0 / 10e12

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15

fig=figure(2, figsize = (6, 8))
clf()
xlabel("k")
ylabel("w/w0")
plot(w/w0,color="black")
ylim((0.0,3))
xlim((0,size(w)[1]))
figure(2, figsize = (6,8))