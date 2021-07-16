import phonopy
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import os
from pymatgen.io.vasp.inputs import Poscar
from ___constants_vasp import VASP_FREQ_TO_INVCM_UNITS
from ___constants_output import DEFAULT_MODE_LABELS
from __directory_searchers import checkPath
from ___helpers_parsing import succ
from random import randint


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
        self.eigensys = [LA.eig(DM) for DM in self.DM_at_Gamma]
        self.evals_real = np.real(np.array([v[0] for v in self.eigensys]))
        self.evecs = np.array([v[1] for v in self.eigensys])

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
            ax.set_title(f"Phonons of mode {modeidx}" + r'at $\Gamma$')
            this_outname = outname[:outname.index('.')] + f'_{modeidx}' + outname[outname.index('.'):]
            fig.colorbar(cf, shrink=0.43, pad=0.05)
            cf.set_clim(min(colors), max(colors))
            fig.savefig(self.outdir + this_outname)
        succ(f"Successfully plotted configuration maps of mode indices {modeidxs}")


class TwistedRealspacePhonon:
    """
    Plots twisted phonons at Gamma point in realspace. Modes are found as u(k=Gamma|G), 
    Fourier transformed into u(k=Gamma|b), then plotted as pairs {(r=A_delta(b), u(k=Gamma|b))}
    in realspace.
    """
    def __init__(self, b_matrix, ph_list, A_delta):
        self.b_matrix = b_matrix[:,:2]
        self.r_matrix = (LA.inv(A_delta) @ b_matrix.T).T
        self.DM_at_Gamma = [ph.get_dynamical_matrix_at_q([0,0,0]) for ph in ph_list]
        self.__diagonalize_DMs()
    def __diagonalize_DMs(self):
        self.eigensys = [np.eig(DM) for DM in self.DM_at_Gamma]
        self.evals = [v[0] for v in self.eigensys]
        self.evecs = [v[1] for v in self.eigensys]
    def __phonon_inverse_fourier(self):
        pass
    def plot_phonons(self):
        pass