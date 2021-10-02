import numpy as np
import hiphive
from hiphive import ForceConstants, ClusterSpace, ForceConstantPotential
from hiphive import enforce_rotational_sum_rules
from hiphive.utilities import extract_parameters
from matplotlib import pyplot as plt
from ase.io import read
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

THz_to_meV = 4.13567

import phonopy
def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])

def plot_dispersion(fcs, color, label):
    # set force constants
    phon.set_force_constants(fcs.get_fc_array(order=2))

    # get dispersion
    phon.set_band_structure(bands)
    _, qnorms, freqs, _ = phon.get_band_structure()
    qnorms = np.hstack(qnorms)
    freqs = THz_to_meV * np.vstack(freqs)
    lines = ax1.plot(qnorms, freqs, color=color)
    lines[0].set_label(label)

ROOT = './'
# define band path
N_band_pts = 500
G2K = get_band(np.array([0.0, 0.0, 0.0]), np.array([1/3, 1/3, 0.0]), N_band_pts)
bands = [G2K]

# monolayer Janus
SP_PATH = ROOT+"SPOSCAR"
PU_PATH = ROOT+"POSCAR_unit"
FC_PATH = ROOT+"FORCE_CONSTANTS"

ph = phonopy.load(supercell_filename=SP_PATH,
                  force_constants_filename=FC_PATH)
supercell = read(SP_PATH)
prim = read(PU_PATH)
fcs_phonopy = ForceConstants.read_phonopy(supercell, FC_PATH)

# parameters
cutoff = 4
dim = 6 # Looks like this has to match with the supercell dimensions

# map fcs_phonopy onto ClusterSpace and produce fcs_hiphive
cs = ClusterSpace(prim, [cutoff])
parameters = extract_parameters(fcs_phonopy, cs)
fcp = ForceConstantPotential(cs, parameters)
fcs_hiphive = fcp.get_force_constants(supercell)

# enforce rotational sum rules and produce fcs_hiphive_rot
enforced_parameters = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang'])
fcp_rot = ForceConstantPotential(cs, enforced_parameters)
fcs_hiphive_rot = fcp_rot.get_force_constants(supercell)

phonopy_prim = PhonopyAtoms(numbers=prim.numbers, positions=prim.positions, cell=prim.cell)
phon = Phonopy(phonopy_prim, supercell_matrix=np.diag([dim, dim, 1]), primitive_matrix=None)

# plotting
fig = plt.figure(figsize=(5.5, 3.8))
ax1 = fig.add_subplot(111)

ax1.axhline(y=0.0, ls='-', c='k', lw=1.0)
plot_dispersion(fcs_phonopy, 'tab:blue', 'phonopy raw')
plot_dispersion(fcs_hiphive, 'tab:orange', 'phonopy-based FCP')
plot_dispersion(fcs_hiphive_rot, 'tab:green', 'phonopy-based FCP rotational invariant')

# show
ax1.legend(loc=3)
ax1.set_xlim(0, 0.1)
ax1.set_ylim(-8, 15)

ax1.set_xlabel('Wave vector')
ax1.set_ylabel('Energy (meV)')

fig.tight_layout()
fig.savefig(ROOT+'hiphive.pdf')
