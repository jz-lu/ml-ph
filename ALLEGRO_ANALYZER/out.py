import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from __class_DataOutput import DataOutput
ABS_MIN_ENERGY = -40.02975439 # energy of minimum configuration to shift everything against.

a = []
a.append((np.asarray([0., 0., 0.]), 0.226122, -36.9177224))
a.append((np.asarray([0.      , 0.333333, 0.      ]), 0.22206, -36.9179183))
a.append((np.asarray([0.      , 0.666667, 0.      ]), 0.221466, -36.91792992))
a.append((np.asarray([0.333333, 0.      , 0.      ]), 0.221393, -36.91795228))
a.append((np.asarray([0.333333, 0.333333, 0.      ]), 0.21728, -36.9179966))
a.append((np.asarray([0.333333, 0.666667, 0.      ]), 0.223795, -36.91788055))
a.append((np.asarray([0.666667, 0.      , 0.      ]), 0.223927, -36.91789084))
a.append((np.asarray([0.666667, 0.333333, 0.      ]), 0.223875, -36.91787807))
a.append((np.asarray([0.666667, 0.666667, 0.      ]), 0.21511, -36.91800157))

p = Poscar.from_file('POS_TEST')
lat = p.as_dict()['structure']['lattice']['matrix'][:-1]
cob = np.transpose([i[:-1] for i in lat])
print("\n\nCOB:", cob, "\n\n")

d = DataOutput('./out', a, cob, ABS_MIN_ENERGY)
d.output_all_analysis()