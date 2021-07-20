import phonopy
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import os
from itertools import product as prod
from random import randint
from pymatgen.io.vasp.inputs import Poscar
from ___constants_vasp import VASP_FREQ_TO_INVCM_UNITS
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
        assert 0 <= shift < self.dms[0], f"Invalid shift index {shift}, max={self.dms[0]}"
        assert len(modeidxs) == len(labels), f"Labels and mode indices must be in 1-1 correspondence, but is {labels} and {modeidxs}"
        for modeidx, lab in zip(modeidxs, labels):
            plt.clf(); fig = plt.figure()
            ax = fig.gca(projection='3d')
            this_outname = outname[:outname.index('.')] + f'_{shift}_{modeidx}' + outname[outname.index('.'):]
            wf = np.real(self.evecs[shift, :, modeidx])
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


"""
Plots twisted phonons at Gamma point in realspace under continuum model. Since the twisted DM
is expressed in Fourier basis, phonons are found as u(k|GM), Fourier transformed into 
u(k=Gamma|r), then plotted as pairs {(r, u(k|r))} in realspace, with r sampled in a supercell.

Note: while this was intended for analyzing phonons at the Gamma point, it makes 
no assumption of what k is, and thus can be used for arbitrary k (just pass into `DM_at_Gamma`).
"""
class TwistedRealspacePhonon:
    def __init__(self, theta, GM_set, DM_at_Gamma, n_at, poscar_uc : Poscar, gridsz=11, outdir='.', cut=6):
        self.GM_set = GM_set; self.n_G = len(GM_set)
        self.cut = cut; self.n_at = n_at
        self.DM_at_Gamma = DM_at_Gamma
        angle = np.deg2rad(theta)
        A_delta2inv = LA.inv(np.array([[1-np.cos(angle), -np.sin(angle)],[np.sin(angle), 1-np.cos(angle)]]))
        self.sc_lattice = A_delta2inv @ poscar_uc.structure.lattice.matrix[:2,:2]
        x = np.linspace(0, 1, num=gridsz, endpoint=False)
        self.d = 3; self.theta = theta
        self.gridsz = gridsz; self.n_r = gridsz**2
        self.r_matrix = np.array(list(prod(x, x))); assert self.r_matrix.shape == (self.n_r, 2)
        self.r_matrix = (self.sc_lattice @ self.r_matrix.T).T # make Cartesian
        self.outdir = checkPath(os.path.abspath(outdir))
        self.phtnsr = None; self.rphtnsr = None # Fourier and real space phonons
        self.rphtnsr_shape = (self.n_r, self.n_at, self.cut, self.d) # realspace phonons
        self.phtnsr_shape = (self.n_G, self.n_at, self.cut, self.d)
        print("Initializing twisted realspace phonon object...")
        self.__phonon_inverse_fourier()
        print(f"Twisted realspace phonon object initialized. Cut={cut}")

    def __diagonalize_DMs(self):
        print("Diagonalizing and sorting dynamical matrix at Gamma...")
        def sorted_eigsys(A):
            vals, vecs = LA.eig(A); idxs = vals.argsort()   
            vals = vals[idxs]; vecs = vecs[:,idxs]
            return vals, vecs
        evals, evecs = sorted_eigsys(self.DM_at_Gamma)
        print("Diagonalized.")
        print("Transforming eigenmatrix into truncated Fourier phonon tensor...")
        # In this case n_at refers to number of atoms in entire bilayer unit cell
        evecs = np.array(evecs[:,:self.cut]).T # shape: (C, n_G x n_at x d) [c := cut = num evecs]
        evecs = np.array(np.split(evecs, self.n_G, axis=1)) # shape: (n_G, C, n_at x d)
        evecs = np.array(np.split(evecs, self.n_at, axis=2)) # shape: (n_at, n_G, C, d)
        self.phtnsr = np.transpose(evecs, axes=(1,0,2,3)) # shape: (n_G, n_at, C, d)
        assert self.phtnsr.shape == self.phtnsr_shape, f"Unexpected phonon tensor shape {self.phtnsr.shape}, expected {self.phtnsr_shape}"
        self.evals_real = np.real(np.array(evals)[:self.cut])
        print("Fourier phonon tensor constructed.")

    def __build_modes(self):
        self.__diagonalize_DMs()
        print("Signing and unitizing phonon modes...")
        self.modes = np.sign(self.evals_real) * np.sqrt(np.abs(self.evals_real)) * VASP_FREQ_TO_INVCM_UNITS
        print("Modes prepared.")

    def __phonon_inverse_fourier(self):
        self.__build_modes()
        print("Building realspace phonon tensor...")
        self.rphtnsr = np.array([sum([G_blk * np.exp(-1j * np.dot(GM, r)) for G_blk, GM in zip(self.phtnsr, self.GM_set)]) for r in self.r_matrix])
        assert self.rphtnsr.shape == self.rphtnsr_shape, f"Unexpected phonon tensor shape {self.rphtnsr.shape}, expected {self.rphtnsr_shape}"
        self.rphtnsr = np.transpose(self.rphtnsr, axes=(1,2,0,3)) # shape: (n_at, C, n_r, d)
        assert self.rphtnsr.shape == (self.n_at, self.cut, self.n_r, self.d)
        print("Realspace phonon tensor built.")

    def plot_phonons(self, outname='phreal.png'):
        coords = self.r_matrix
        np.save(self.outdir + "rphtnsr.npy", self.rphtnsr)
        for at_i, atomic_blk in enumerate(self.rphtnsr):
            for m_j, phonons in enumerate(atomic_blk):
                # phonons = np.array([wf/LA.norm(wf) for wf in phonons])
                plt.clf(); fig = plt.figure(); ax = fig.gca(projection='3d')
                plt.quiver(coords[:,0], coords[:,1], np.zeros(self.n_r), 
                            phonons[:,0], phonons[:,1], phonons[:,2], length=0.5)
                ax.scatter(coords[:,0], coords[:,1], np.zeros(self.n_r), c='black')
                this_outname = outname[:outname.index('.')] + f'_{m_j}_{at_i}' + outname[outname.index('.'):]
                plt.title(f"Normalized phonons (" + r"$\theta=$" + '%.1lf'%self.theta + r"$^\circ$" + f", mode {m_j}, atom {at_i}) in supercell")
                fig.savefig(self.outdir + this_outname)
                plt.close(fig)
                update(f"Wrote twisted phonons in realspace to {self.outdir + this_outname}")
        succ(f"Successfully generated {self.cut * self.n_at} realspace twisted phonon plots")

