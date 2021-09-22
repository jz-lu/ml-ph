import numpy as np
from matplotlib import pyplot as plt
from ase.io import read
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from ase.lattice.hexagonal import Graphene

from hiphive import ForceConstants, ClusterSpace, ForceConstantPotential
from hiphive import enforce_rotational_sum_rules
from hiphive.utilities import extract_parameters


THz_to_meV = 4.13567


def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])


def plot_dispersion(fcs, color, label):
    # set force constants
    phonopy.set_force_constants(fcs.get_fc_array(order=2))

    # get dispersion
    phonopy.set_band_structure(bands)
    _, qnorms, freqs, _ = phonopy.get_band_structure()
    qnorms = np.hstack(qnorms)
    freqs = THz_to_meV * np.vstack(freqs)
    lines = ax1.plot(qnorms, freqs, color=color)
    lines[0].set_label(label)


# parameters
cutoff = 8
dim = 7
mesh = [24] * 3
a0 = 2.466340583
vacuum = 40

# set up prim structure
prim = Graphene(symbol='C', latticeconstant={'a': a0, 'c': vacuum})

# define band path
N_band_pts = 500
G2K = get_band(np.array([0.0, 0.0, 0.0]), np.array([1/3, 1/3, 0.0]), N_band_pts)
bands = [G2K]

# init phonopy
phonopy_prim = PhonopyAtoms(numbers=prim.numbers, positions=prim.positions, cell=prim.cell)
phonopy = Phonopy(phonopy_prim, supercell_matrix=np.diag([dim, dim, 1]), primitive_matrix=None)
phonopy.generate_displacements()

# load the precalculated single phonopy supercell with displacements and forces
supercell = read('graphene_phonopy_supercell.extxyz')

# Add forces to phonopy and calculate force constants
phonopy.set_forces([supercell.arrays['forces']])
phonopy.produce_force_constants()
fcs_phonopy = ForceConstants.from_arrays(supercell, phonopy.get_force_constants())

# map fcs_phonopy onto ClusterSpace and produce fcs_hiphive
cs = ClusterSpace(prim, [cutoff])
parameters = extract_parameters(fcs_phonopy, cs)
fcp = ForceConstantPotential(cs, parameters)
fcs_hiphive = fcp.get_force_constants(supercell)

# enforce rotational sum rules and produce fcs_hiphive_rot
enforced_parameters = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang'])
fcp_rot = ForceConstantPotential(cs, enforced_parameters)
fcs_hiphive_rot = fcp_rot.get_force_constants(supercell)


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
fig.savefig('graphene_phonon_dispersion.pdf')