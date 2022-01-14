import phonopy
import numpy as np
import numpy.linalg as LA
from scipy import linalg as SLA
from math import floor, log10, sqrt
import matplotlib.pyplot as plt
import os
from itertools import product as prod
from random import randint
from pymatgen.io.vasp.inputs import Poscar
from ___constants_vasp import VASP_FREQ_TO_INVCM_UNITS, Z_LAYER_SEP
from ___constants_output import DEFAULT_MODE_LABELS
from __directory_searchers import checkPath
from ___helpers_parsing import succ, update

"""
Phonon modes in configuration space; phonon displacements for each atom in realspace.
"""
class PhononConfig:
    def __init__(self, b_matrix, cob, ph_list, outdir='.'):
        self.outdir = checkPath(os.path.abspath(outdir))
        self.DM_at_Gamma = [ph.get_dynamical_matrix_at_q([0,0,0]) for ph in ph_list]
        self.dms = dms = self.DM_at_Gamma[0].shape
        self.ncfg = len(ph_list)
        assert len(dms) == 2 and dms[0] == dms[1] and dms[0] % 3 == 0, f"Invalid force constants matrix shape {dms}"
        print(f"Bilayer dynamical matrix in configuration space has shape {dms}")
        self.b_matrix = (cob @ b_matrix[:,:2].T).T
        self.__build_modes()
        print("Initialized PhononConfig object to analyze modes at Gamma point")

    def __diagonalize_DMs(self):
        def sorted_eigsys(A):
            vals, vecs = LA.eig(A); idxs = vals.argsort()   
            vals = vals[idxs]; vecs = vecs[:,idxs]
            return vals, vecs
        eigsys = [sorted_eigsys(DM) for DM in self.DM_at_Gamma]
        evals = [v for v,_ in eigsys]; evecs = [v for _,v in eigsys]
        self.evecs = np.array(evecs)
        self.evals_real = np.real(np.array(evals))

    def __build_modes(self):
        self.__diagonalize_DMs()
        self.modes = np.sign(self.evals_real) * np.sqrt(np.abs(self.evals_real)) * VASP_FREQ_TO_INVCM_UNITS

    def plot_mode_quiver(self, poscar : Poscar, shift=0, modeidxs=np.arange(6), labels=DEFAULT_MODE_LABELS, outname='cfgquiver.png'):
        coords = poscar.structure.cart_coords
        assert 0 <= shift < self.ncfg, f"Invalid shift index {shift}, max={self.ncfg}"
        assert len(modeidxs) == len(labels), f"Labels and mode indices must be in 1-1 correspondence, but is {labels} and {modeidxs}"
        for modeidx, lab in zip(modeidxs, labels):
            plt.clf(); fig = plt.figure()
            ax = fig.gca(projection='3d')
            this_outname = outname[:outname.index('.')] + f'_{shift}_{modeidx}' + outname[outname.index('.'):]
            wf = np.real(self.evecs[shift, :, modeidx])
            print(f"WF-{shift}-{modeidx}: {wf}")
            ax.scatter(coords[:,0], coords[:,1], coords[:,2])
            plt.quiver(coords[:,0], coords[:,1], coords[:,2], 
                       wf[0::3], wf[1::3], wf[2::3], length=1, normalize=True)
            plt.title(f"Phonons (mode {modeidx}) at shift {shift}")
            fig.savefig(self.outdir + this_outname)
        succ(f"Successfully plotted quivers of mode indices {modeidxs} at shift index {shift}")
            
    def plot_mode_config(self, modeidxs=np.arange(6), labels=DEFAULT_MODE_LABELS, outname='cfgph.png'):
        for modeidx, lab in zip(modeidxs, labels):
            plt.clf(); fig, ax = plt.subplots()
            colors = self.modes[:,modeidx]
            cf = ax.scatter(self.b_matrix[:,0], self.b_matrix[:,1], s=300, c=colors, cmap='RdYlBu')
            ax.set_aspect('equal')
            ax.set_title(f"Phonons of mode {modeidx} " + r'at $\Gamma$')
            this_outname = outname[:outname.index('.')] + f'_{modeidx}' + outname[outname.index('.'):]
            fig.colorbar(cf, shrink=0.43, pad=0.05)
            cf.set_clim(min(colors), max(colors))
            fig.savefig(self.outdir + this_outname)
        succ(f"Successfully plotted configuration maps of mode indices {modeidxs}")

def make_unit_square():
    """Helper function that returns the coordinates of a unit square"""
    npts = 11
    a = np.linspace(0, 1, npts, endpoint=True)
    a0 = np.zeros_like(a); a1 = np.ones_like(a)
    l1 = np.stack((a, a0), axis=1); l2 = np.stack((a1, a), axis=1)
    l3 = np.stack((np.flip(a), a1), axis=1); l4 = np.stack((a0, np.flip(a)), axis=1)
    sq = np.concatenate((l1, l2, l3, l4)); assert sq.shape == (4*npts, 2)
    return sq

"""
Plots twisted phonons at Gamma point in realspace under continuum model. Since the twisted DM
is expressed in Fourier basis, phonons are found as u(k|GM), Fourier transformed into 
u(k=Gamma|r), then plotted as pairs {(r, u(k|r))} in realspace, with r sampled in a supercell.

Note: while this was intended for analyzing phonons at the Gamma point, it makes 
no assumption of what k is, and thus can be used for arbitrary k (just pass into `DM_at_Gamma`).

Bugs: fails if number of atoms in each layer differ. In practice this is not much of a concern
anyway, since the lattice constants will be too far apart to form a stable bilayer. If for
some reason it is useful for different numbers of atoms, simply adjust layer slicing in formation 
of `phtnsr` to slice proportional to ratio of atoms, instead of in half.
"""
class TwistedRealspacePhonon:
    def __init__(self, theta, k, GM_set, DMset_at_k, n_at, bl_masses, 
                 poscars_uc, gridsz=13, outdir='.', modeidxs=np.linspace(0,6,6), kpt=r'$\Gamma$',
                 RSPC_SUPERCELL_SIZE=2):
        assert len(bl_masses) == n_at
        self.RSPC_SUPERCELL_SIZE = RSPC_SUPERCELL_SIZE
        self.kpt = kpt
        self.bl_masses = bl_masses
        self.GM_set = GM_set; self.n_G = len(GM_set); self.k = k
        self.modeidxs = np.array(modeidxs).astype(int)
        self.nmodes = len(modeidxs); assert self.nmodes >= 1
        self.n_at = n_at
        self.DMset_at_k = DMset_at_k
        angle = np.deg2rad(theta)
        A_delta2inv = LA.inv(np.array([[1-np.cos(angle), np.sin(angle)],[-np.sin(angle), 1-np.cos(angle)]]))
        A0 = poscars_uc[0].structure.lattice.matrix[:2,:2].T
        self.sc_lattice = A_delta2inv @ A0
        self.at_pos = np.concatenate([p.structure.cart_coords for p in poscars_uc], axis=0)
        self.at_pos[self.n_at//2:,2] += Z_LAYER_SEP*poscars_uc[0].structure.lattice.matrix[-1,-1] # make interlayer space
        x = np.linspace(-RSPC_SUPERCELL_SIZE/2, RSPC_SUPERCELL_SIZE/2, num=gridsz*RSPC_SUPERCELL_SIZE, endpoint=False)
        self.d = 3; self.theta = theta
        self.gridsz = gridsz; self.n_r = (gridsz * RSPC_SUPERCELL_SIZE)**2
        self.r_matrix = np.array(list(prod(x, x))); assert self.r_matrix.shape == (self.n_r, 2)
        self.diagidxs = self.r_matrix[:,0] == self.r_matrix[:,1]
        self.r_linecut_direct = self.r_matrix[self.diagidxs]
        self.r_matrix = self.r_matrix @ self.sc_lattice.T # make Cartesian
        self.moire_boundary = make_unit_square() @ self.sc_lattice.T # boundary line of moire cell
        self.outdir = checkPath(os.path.abspath(outdir))
        self.phtnsr = None; self.rphtnsr = None # Fourier and real space phonons
        self.old_rphtnsr_shape = (self.n_r, self.n_at, self.nmodes, self.d) # realspace phonons
        self.phtnsr_shape = (self.n_G, self.n_at, self.nmodes, self.d)
        print("Initializing twisted realspace phonon object...")
        self.__phonon_inverse_fourier()
        print(f"Twisted realspace phonon object initialized. Mode indices={self.modeidxs}")

    def __DM_to_phtnsr(self):
        print("Diagonalizing and sorting dynamical matrix at Gamma...")
        def extract_G0_eig(A):
            vals, vecs = LA.eig(A)
            vecs_by_layer = np.split(vecs, 2, axis=0)
            # breakpoint()
            for li in range(2):
                vecs_by_layer[li] = np.array(np.split(vecs_by_layer[li], self.n_G, axis=0)[0]) # G0 only
            # breakpoint()
            return (vals, vecs_by_layer)
        def sorted_filtered_eigsys(M_set):
            assert len(M_set) == self.n_G
            raw_eigsys = [extract_G0_eig(M) for M in M_set]
            vals = np.real(raw_eigsys[0][0]) # the eigenvalues corr to bands are at k = k + G0
            
            # In this case n_at refers to number of atoms in entire bilayer unit cell
            vecs = np.array([vec for _,vec in raw_eigsys]) # shape: (n_G, 2, n_at/2 x d, num_evecs=n_G x n_at x d)
            idxs = vals.argsort()  
            vals = vals[idxs]; vecs = vecs[:,:,:,idxs]
            vecs = np.array(np.split(vecs, self.n_at//2, axis=2)) # shape: (n_at/2, n_G, 2, d, num_evecs)
            vecs = np.transpose(vecs, axes=(2,0,1,4,3)) # shape: (2, n_at/2, n_G, num_evecs, d)
            vecs = np.concatenate(vecs, axis=0) # shape: (n_at, n_G, num_evecs, d)
            assert vals.shape[0] == vecs.shape[2]
            return vals, vecs

        evals, evecs = sorted_filtered_eigsys(self.DMset_at_k)
        self.evals_real = evals[self.modeidxs]
        evecs = evecs[:,:,self.modeidxs] # shape: (n_at, n_G, C, d) [C := num evecs user wants]
        print("Diagonalized.")
        print("Transforming eigenmatrix into truncated Fourier phonon tensor...")
        
        assert evecs.shape == (self.n_at, self.n_G, self.nmodes, self.d), f"Incorrect shape {evecs.shape}"
        self.phtnsr = np.transpose(evecs, axes=(1,0,2,3)) # shape: (n_G, n_at, C, d)
        assert self.phtnsr.shape == self.phtnsr_shape, f"Unexpected phonon tensor shape {self.phtnsr.shape}, expected {self.phtnsr_shape}"
        print(f"Fourier phonon tensor constructed: shape {self.phtnsr.shape}")

    def __build_modes(self):
        self.__DM_to_phtnsr()
        print("Signing and unitizing phonon modes...")
        self.modes = np.sign(self.evals_real) * np.sqrt(np.abs(self.evals_real)) * VASP_FREQ_TO_INVCM_UNITS
        print("Modes prepared.")

    def __phonon_inverse_fourier(self):
        self.__build_modes()
        print(f"Building realspace phonon tensor...(k={self.k}, |k| = {LA.norm(self.k)})")
        
        self.rphtnsr = np.array([sum([G_blk * np.exp(-1j * np.dot(GM + self.k, r)) 
            for G_blk, GM in zip(self.phtnsr, self.GM_set)]) for r in self.r_matrix])
        assert self.rphtnsr.shape == self.old_rphtnsr_shape, \
            f"Unexpected phonon tensor shape {self.rphtnsr.shape}, expected {self.old_rphtnsr_shape}"
            
        self.rphtnsr = np.transpose(self.rphtnsr, axes=(1,2,0,3)) # shape: (n_at, C, n_r, d)
        assert self.rphtnsr.shape == (self.n_at, self.nmodes, self.n_r, self.d)
        print(f"Realspace phonon tensor built: shape {self.rphtnsr.shape}")
        self.mnormed_tnsr = np.array([y/sqrt(x) for x, y in zip(self.bl_masses, self.rphtnsr)])
        # breakpoint()
    
    def plot_spatial_avgs(self, outname='spavg.png'):
        avgtnsr = np.mean(self.rphtnsr, axis=2) # shape: (n_at, C, d)
        avgtnsr = np.transpose(avgtnsr, axes=(1,0,2)) # shape: (C, n_at, d)
        np.save(self.outdir + f"avgtnsr_{self.kpt}.npy", avgtnsr)
        print(f"Realspace average tensor shape: {avgtnsr.shape}")
        nuc_at = self.n_at//2
        plt.clf(); plt.plot(avgtnsr[0,:,0], avgtnsr[0,:,1]) # just to get the right dom/range
        left, right = plt.xlim()
        x = np.zeros(self.n_at); y = np.linspace(10*left,10*right,num=self.n_at); yset = False; lim = None
        lcol = ['black']*nuc_at + ['maroon']*nuc_at
        for m_j, phonons in enumerate(avgtnsr): # shape: (n_at, d)
            plt.clf(); fig, axes = plt.subplots(nrows=1, ncols=2)
            if not yset:
                dummy = axes[0].quiver(x, y, phonons[:,0], phonons[:,1], color=lcol)
                lim = np.array(axes[0].get_xlim()); lim = lim/2
                y = np.linspace(lim[0]*2/3, lim[1]*2/3, num=self.n_at)
                dummy.remove()
            axes[0].quiver(x, y, phonons[:,0], phonons[:,1], color=lcol) # xy projection
            lscale = floor(log10(np.mean(list(map(LA.norm, phonons[:,0:2])))))
            axes[0].text(lim[0], lim[1], str(lscale))
            axes[1].quiver(x, y, phonons[:,0], phonons[:,2], color=lcol) # xz projection
            lscale = floor(log10(np.mean(list(map(LA.norm, phonons[:,1:3])))))
            axes[0].text(lim[0], lim[1], str(lscale))
            axes[0].set_title('xy')
            axes[1].set_title('xz')
            for ax in axes:
                ax.set_ylim(lim)
                ax.set_xlim(lim)
                ax.set_aspect('equal')
            plt.suptitle(f"Mean field phonons from continuum (mode {self.modeidxs[m_j]})")
            this_outname = outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_k{self.kpt}' + outname[outname.index('.'):]
            fig.savefig(self.outdir + this_outname)
            plt.close(fig)
        succ("Successfully outputted spatial-average phonon plots")
    
    def plot_atomic_avgs(self, outname='atavg.png'):
        # Take a diagonal cut of the sampling, and plot averages over all atoms
        # Key idea: shear/LB modes will cancel on average, while translations will not.
        # For this function only n_r := number of realsoace points along the line cut.
        avgtnsr = np.mean(self.rphtnsr, axis=0) # shape: (C, n_r, d)
        avgtnsr = avgtnsr[:,self.diagidxs,:] # keep diagonal terms only
        print(f"Atomic average tensor shape: {avgtnsr.shape}")
        x = self.r_linecut_direct[:,0] / 15; y = self.r_linecut_direct[:,1] / 15
        assert len(x) == sum(self.diagidxs)
        for m_j, phonons in enumerate(avgtnsr): # shape: (n_r, d)
            plt.clf(); fig, axes = plt.subplots(nrows=1, ncols=2)
            axes[0].quiver(x, y, phonons[:,0], phonons[:,1]) # xy projection
            axes[1].quiver(x, y, phonons[:,0], phonons[:,2]) # xz projection
            axes[0].set_title('xy')
            axes[1].set_title('xz')
            for ax in axes:
                ax.set_aspect('equal')
            plt.suptitle(f"Mean atom phonons in diagonal cut (mode {self.modeidxs[m_j]})")
            this_outname = outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_k{self.kpt}' + outname[outname.index('.'):]
            fig.savefig(self.outdir + this_outname)
            plt.close(fig)
        succ("Successfully outputted atomic-average phonon plots")

    # Generate a plot of the phonons averaged over each layer;
    # there is one plot per mode, per layer
    def plot_phonons(self, outname='phreal.png', zcolmesh=False):
        coords = self.r_matrix
        np.save(self.outdir + f"rphtnsr_k{self.kpt}.npy", self.rphtnsr)

        layer_blks = np.split(self.mnormed_tnsr, 2, axis=0) # split back by layer, then avg it
        layer_blks = np.abs(list(map(lambda x: np.mean(x, axis=0), layer_blks)))
        assert layer_blks.shape == (2, self.nmodes, self.n_r, self.d)
        zbound = np.max(np.abs(layer_blks[:,:,:,2]))
        for l_i, layer_blk in enumerate(layer_blks):
            l_i += 1 # index layers by 1
            for m_j, phonons in enumerate(layer_blk):
                phonons = np.real(phonons) # just take the real component
                z = phonons[:,2]
                plt.clf(); fig, ax = plt.subplots(figsize=(3.5*self.RSPC_SUPERCELL_SIZE, 5.5*self.RSPC_SUPERCELL_SIZE))
                plt.rc('font', size=8*self.RSPC_SUPERCELL_SIZE)
                ax.plot(self.moire_boundary[:,0], self.moire_boundary[:,1], c="limegreen", alpha=0.8)
                plt.quiver(coords[:,0], coords[:,1],    # positions
                            phonons[:,0], phonons[:,1], # arrows
                            z,                          # arrow colors
                            cmap='CMRmap')
                (xm, xp), (ym, yp) = plt.xlim(), plt.ylim()
                max_xy = np.max([LA.norm(phonon[:-1]) for phonon in phonons])
                max_z = np.max(np.abs(z))
                ax.text(0.02*(xp-xm)+xm, 0.02*(yp-ym)+ym, r'$\delta u_{xy} = %.3E$'%max_xy)
                ax.text(0.06*(xp-xm)+xm, 0.06*(yp-ym)+ym, r'$\delta u_{z} = %.3E$'%max_z)
                ax.text(0.10*(xp-xm)+xm, 0.10*(yp-ym)+ym, r'$\omega = %.3f$'%self.modes[m_j])
                plt.xlabel("x"); plt.ylabel("y")
                ax.scatter(coords[:,0], coords[:,1], c='black', s=0.2)
                ax.set_aspect('equal')
                fname = self.kpt[2:-1] if self.kpt[0] == "$" else self.kpt
                this_outname = outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_{l_i}_k-{fname}' + outname[outname.index('.'):]
                plt.title(r"$\theta=$" + '%.1lf'%self.theta + r"$^\circ,$" + f" Mode {self.modeidxs[m_j]}, Layer {l_i} at " + self.kpt)
                plt.colorbar(shrink=0.5)
                plt.clim(-zbound, zbound)
                fig.savefig(self.outdir + this_outname)
                plt.close(fig)

        if zcolmesh:
            for l_i, layer_blk in enumerate(layer_blks):
                l_i += 1 # index layers by 1
                for m_j, phonons in enumerate(layer_blk):
                    phonons = np.real(phonons) # just take the real component
                    z = phonons[:,2]
                    plt.clf(); fig, ax = plt.subplots(figsize=(3.5*self.RSPC_SUPERCELL_SIZE, 5.5*self.RSPC_SUPERCELL_SIZE))
                    plt.rc('font', size=8*self.RSPC_SUPERCELL_SIZE)
                    plt.tricontourf(coords[:,0], coords[:,1], z, cmap='CMRmap', levels=201)
                    ax.plot(self.moire_boundary[:,0], self.moire_boundary[:,1], c="limegreen", alpha=0.8)
                    (xm, xp), (ym, yp) = plt.xlim(), plt.ylim()
                    max_xy = np.max([LA.norm(phonon[:-1]) for phonon in phonons])
                    max_z = np.max(np.abs(z))
                    ax.text(0.02*(xp-xm)+xm, 0.02*(yp-ym)+ym, r'$\delta u_{xy} = %.3E$'%max_xy)
                    ax.text(0.06*(xp-xm)+xm, 0.06*(yp-ym)+ym, r'$\delta u_{z} = %.3E$'%max_z)
                    ax.text(0.10*(xp-xm)+xm, 0.10*(yp-ym)+ym, r'$\omega = %.3f$'%self.modes[m_j])
                    plt.xlabel("x"); plt.ylabel("y")
                    ax.set_aspect('equal')
                    fname = self.kpt[2:-1] if self.kpt[0] == "$" else self.kpt
                    this_outname = "COL_" + outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_{l_i}_k-{fname}' + outname[outname.index('.'):]
                    plt.title(r"$\theta=$" + '%.1lf'%self.theta + r"$^\circ,$" + f" Mode {self.modeidxs[m_j]}, Layer {l_i} at " + self.kpt)
                    cb = plt.colorbar(shrink=0.5)
                    cb.mappable.set_clim(-zbound, zbound)
                    fig.savefig(self.outdir + this_outname)
                    plt.close(fig)

                update(f"Wrote twisted phonons in realspace to {self.outdir + this_outname}")
        succ(f"Successfully generated {self.nmodes * self.n_at} realspace twisted phonon plots")
        
    # One plot per mode, per atom, per layer
    def plot_phonons_per_atom(self, outname='phat.png', zcolmesh=False):
        coords = self.r_matrix
        np.save(self.outdir + f"rphtnsr_k{self.kpt}.npy", self.rphtnsr)

        mnormed_tnsr = np.array([y/sqrt(x) for x, y in zip(self.bl_masses, self.rphtnsr)])
        layer_blks = np.array(np.split(mnormed_tnsr, 2, axis=0)) # split back by layer, then avg it
        zbound = np.max(np.abs(layer_blks[:,:,:,:,2]))
        for l_i, layer_blk in enumerate(layer_blks):
            l_i += 1 # index layers by 1
            for at_k, at_blk in enumerate(layer_blk):
                for m_j, phonons in enumerate(at_blk):
                    phonons = np.real(phonons) # just take the real component
                    z = phonons[:,2]
                    plt.clf(); fig, ax = plt.subplots(figsize=(3.5*self.RSPC_SUPERCELL_SIZE, 5.5*self.RSPC_SUPERCELL_SIZE))
                    plt.rc('font', size=8*self.RSPC_SUPERCELL_SIZE)
                    ax.plot(self.moire_boundary[:,0], self.moire_boundary[:,1], c="limegreen", alpha=0.8)
                    plt.quiver(coords[:,0], coords[:,1],    # positions
                                phonons[:,0], phonons[:,1], # arrows
                                z,                          # arrow colors
                                cmap='CMRmap')
                    (xm, xp), (ym, yp) = plt.xlim(), plt.ylim()
                    max_xy = np.max([LA.norm(phonon[:-1]) for phonon in phonons])
                    max_z = np.max(np.abs(z))
                    ax.text(0.02*(xp-xm)+xm, 0.02*(yp-ym)+ym, r'$\delta u_{xy} = %.3E$'%max_xy)
                    ax.text(0.06*(xp-xm)+xm, 0.06*(yp-ym)+ym, r'$\delta u_{z} = %.3E$'%max_z)
                    ax.text(0.10*(xp-xm)+xm, 0.10*(yp-ym)+ym, r'$\omega = %.3f$'%self.modes[m_j])
                    plt.xlabel("x"); plt.ylabel("y")
                    ax.scatter(coords[:,0], coords[:,1], c='black', s=0.2)
                    ax.set_aspect('equal')
                    fname = self.kpt[2:-1] if self.kpt[0] == "$" else self.kpt
                    this_outname = outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_{l_i}_{at_k}_k-{fname}' + outname[outname.index('.'):]
                    plt.title(r"$\theta=$" + '%.1lf'%self.theta + r"$^\circ,$" + f" Atom {at_k}, Mode {self.modeidxs[m_j]}, Layer {l_i} at " + self.kpt)
                    plt.colorbar(shrink=0.5)
                    plt.clim(-zbound, zbound)
                    fig.savefig(self.outdir + this_outname)
                    plt.close(fig)

            if zcolmesh:
                for l_i, layer_blk in enumerate(layer_blks):
                    l_i += 1 # index layers by 1
                    for at_k, at_blk in enumerate(layer_blk):
                        for m_j, phonons in enumerate(at_blk):
                            phonons = np.real(phonons) # just take the real component
                            z = phonons[:,2]
                            plt.clf(); fig, ax = plt.subplots(figsize=(3.5*self.RSPC_SUPERCELL_SIZE, 5.5*self.RSPC_SUPERCELL_SIZE))
                            plt.rc('font', size=8*self.RSPC_SUPERCELL_SIZE)
                            plt.tricontourf(coords[:,0], coords[:,1], z, cmap='CMRmap', levels=201)
                            ax.plot(self.moire_boundary[:,0], self.moire_boundary[:,1], c="limegreen", alpha=0.8)
                            (xm, xp), (ym, yp) = plt.xlim(), plt.ylim()
                            max_xy = np.max([LA.norm(phonon[:-1]) for phonon in phonons])
                            max_z = np.max(np.abs(z))
                            ax.text(0.02*(xp-xm)+xm, 0.02*(yp-ym)+ym, r'$\delta u_{xy} = %.3E$'%max_xy)
                            ax.text(0.06*(xp-xm)+xm, 0.06*(yp-ym)+ym, r'$\delta u_{z} = %.3E$'%max_z)
                            ax.text(0.10*(xp-xm)+xm, 0.10*(yp-ym)+ym, r'$\omega = %.3f$'%self.modes[m_j])
                            plt.xlabel("x"); plt.ylabel("y")
                            ax.set_aspect('equal')
                            fname = self.kpt[2:-1] if self.kpt[0] == "$" else self.kpt
                            this_outname = "COL_" + outname[:outname.index('.')] + f'_{self.modeidxs[m_j]}_{l_i}_{at_k}_k-{fname}' + outname[outname.index('.'):]
                            plt.title(r"$\theta=$" + '%.1lf'%self.theta + r"$^\circ,$" + f" Atom {at_k}, Mode {self.modeidxs[m_j]}, Layer {l_i} at " + self.kpt)
                            cb = plt.colorbar(shrink=0.5)
                            cb.mappable.set_clim(-zbound, zbound)
                            fig.savefig(self.outdir + this_outname)
                            plt.close(fig)

                update(f"Wrote per-atom twisted phonons in realspace to {self.outdir + this_outname}")
        succ(f"Successfully generated {self.nmodes * self.n_at} realspace twisted phonon plots")

