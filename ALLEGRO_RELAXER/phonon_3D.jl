
## 3D phonon modes
using FFTW

FFTW.set_num_threads(Sys.CPU_THREADS)
include("Bilayers.jl")
# include("GrapheneParameters.jl")
include("Janus_SSe_180_Parameters.jl")

using PyPlot
using PyCall
using Interpolations
using Optim
using LinearAlgebra

close("all")

## Parameters:
θ = deg2rad(4)

N = 51
hN = div(N,2)
blg = Bilayer(l, θ, K, G)
hull = Hull(blg, N, [0.,0.]);

# Cutofff radius in k-space
G_cutoff = 4 # need to be large for small angles
# number of points to sample on the high symmetry line
nq = 50
nr = 50
mode = "line" # mode = "mesh" or "line", "line" to calculate band structures, "mesh" for DOS (not implemented)
wf_idx = 1

## Here are rotated functionals with symmetry - try *_nosym versions for the simplified functional in the paper
f(u) = Energy(u, hull)
g!(storage, u) = Gradient!(storage, u, hull)
@time results = optimize(f, g!, zeros(2, N^2), LBFGS(), Optim.Options(allow_f_increases=true))

# static solution
u0 = reshape(Optim.minimizer(results), (2, N^2))
# u0 =  blg.E * [1/3; 1/3 ]
∇u = ∇(u0, hull)
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
xlim((-20,20))
ylim((-20,20))
ax=gca()
ax.set_aspect(1)
colorbar()

subplot(2,1,2)
pcolormesh(kx,ky,abs.(uG[2,:,:]))
xlim((-20,20))
ylim((-20,20))
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
evals = zeros(Float64,(size(kline,2),size(q_list,2)*3))

# interlayer
# first calculate the Fourier component of the gradient of the relaxed misfit energy
A = zeros(Float64, (N, N))
plan = plan_bfft(A, (1,2) )

# Fourier decomposition of relaxed GSFE
B = hull.G + 2 * blg.iR * u0 # relaxed configuration space position
VG = zeros(ComplexF64, (3, 3, N, N))

z = reshape(GSFE_height(B, blg), (N, N))
dz = reshape(gradient_GSFE_height(B, blg), (2, N, N))
dz_x = plan * reshape(dz[1, :, :], (N, N))
dz_y = plan * reshape(dz[2, :, :], (N, N))

dz_x = dz_x[ind,:]
dz_x = dz_x[:,ind]
dz_y = dz_y[ind,:]
dz_y = dz_y[:,ind]

VG[3, 1, :, :] = dz_x / N^2
VG[3, 2, :, :] = dz_y / N^2
VG[1, 3, :, :] = dz_x / N^2
VG[2, 3, :, :] = dz_y / N^2

for ci = 1:5
    prefac_tmp = inter_ph(blg, ci)
    gsfe = zeros(ComplexF64, (N, N, 3))

    for j = 1:3
        gsfe[:,:,j] = plan * reshape(GSFE_comp(B, blg, ci, j), (N, N))
    end

    gsfe=gsfe[ind,:,:]
    gsfe=gsfe[:,ind,:] # sort Fourier component of GSFE by G vector

    for c = 1:3
        for a = 1:2, b = 1:2
            global VG[a,b,:,:] = VG[a,b,:,:] - 2 * prefac_tmp[a,b,c] * gsfe[:,:,c] / N^2
        end
    end
end


VG = reshape(VG, (3, 3, N^2))
# VG = VG[:,:,id]

Hinter = zeros(ComplexF64,(size(q_list,2)*3, size(q_list,2)*3))
qidx = zeros(Int64, (size(q_list,2), size(q_list,2)))
dG = zeros(Float64, (2, size(q_list,2), size(q_list,2)))

tol = 1e-6
for q1_idx = 1:size(q_list,2)
    for q2_idx = 1:size(q_list,2)

        id1 = q1_idx*3-2;
        id2 = q2_idx*3-2;
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

            Hinter[id1:id1+2,id2:id2+2]=VG[:,:,qid]

        end
    end
end


a0 = l*sqrt(3)
A0 = sqrt(3)/2*a0^2

evecs_all = zeros(ComplexF64,(size(q_list,2)*3, size(q_list,2)*3, size(kline,2)))

for k_idx = 1:size(kline,2)
    Hintra = zeros(ComplexF64,size(q_list,2)*3, size(q_list,2)*3)

    # an equivalent way to create KG
    # hull = Hull(blg, N, kline_idx[:,k_idx])
    # KG = hull.Hess_Elastic
    # KG = reshape(KG, (2, 2, N, N))
    # KG = KG[:,:,ind,:]
    # KG = KG[:,:,:,ind] # sort KG from smallest G vector to the largest
    # KG = reshape(KG,(2,2,N^2))
    # KG = KG[:,:,id] # within the cutoff radius

    KG = zeros(Float64, (3,3,size(q_list,2)))
    for g_idx = 1:size(KG, 3)
        Gx = q_list[1,g_idx] + kline[1, k_idx]
        Gy = q_list[2,g_idx] + kline[2, k_idx]

        KG[:,:,g_idx] = [(λ+2*μ)*Gx^2+μ*Gy^2  (λ+μ)*Gx*Gy 0
                        (λ+μ)*Gx*Gy (λ+2*μ)*Gy^2+μ*Gx^2 0
                        0 0 2*κ*(Gx^2+Gy^2)]

        id1 = g_idx*3-2
        id2 = id1+2
        Hintra[id1:id2,id1:id2]=KG[:,:,g_idx]

    end

    global H = (Hinter+Hintra)/A0 # convert from eV / cell to eV / A^2

    # global H = Hintra/A0

    global F = eigen(H,sortby=nothing) # default sort by (real(λ), real(λ))
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

w = sqrt.(abs.(evals) * 1e20 * 16.02   ./ rho) # in Hz
lambda = λ/A0; # in eV / A^2
lambda = lambda * 1.602e-19 *1e20; # in J / m^2
LM = a0 / θ *1e-10 # moire length in m
w0 = sqrt(lambda/rho) * (2*pi) / LM # in 1/s
# convert to meV
hbar = 4.1357e-15 # in eV.s
w0_meV = w0 * hbar *1e3
w0_THz = w0 / 10e12

fig=figure(2, figsize = (8, 10))
clf()
plot(w/w0*w0_THz,color="black")
ylim((0.0,1))
figure(2, figsize = (8, 10))

# calculate wavefunction at one k point
# sampling real space positions
# band = range(1, length=10)
# r = EM * [xmesh[:]' .- max(xmesh[:]...)/2;ymesh[:]' .- max(ymesh[:]...)/2]/nq * 4
# rx = reshape(r[1,:], (nr, nr))
# ry = reshape(r[2,:], (nr, nr))
#
# uGx = evecs_all[:,1:2:end,:]
# uGy = evecs_all[:,2:2:end,:]
#
# global ubx = zeros(ComplexF64, (nr^2, size(q_list,2)))
# global uby = zeros(ComplexF64, (nr^2, size(q_list,2)))
# j = complex(0,1)
#
# ph = exp(j * sum(dot.(kline[:,wf_idx], r[:,r_idx]),dims=1)[1])
#
# for r_idx = 1:size(r,2)
#     for q_idx = 1:size(q_list,2)
#         ph_here = j * (r[1,r_idx] * q_list[1, q_idx] + r[2, r_idx] * q_list[2, q_idx])
#         ubx[r_idx, :] += uGx[q_idx, :, wf_idx] .* exp(ph_here) / size(q_list,2) .* ph
#         uby[r_idx, :] += uGy[q_idx, :, wf_idx] .* exp(ph_here) / size(q_list,2) .* ph
#     end
# end
#
#
# ubx = reshape(ubx, (nr, nr, size(ubx,2)))
# uby = reshape(uby, size(ubx))

# band_idx = 2
# fig = figure(5,figsize = (8, 6))
# clf()
# quiver(rx,ry,real.(ubx[:,:,band_idx]),real.(uby[:,:,band_idx]))
# ax=gca()
# ax.set_aspect(1)
# figure(5)
