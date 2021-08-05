from pymatgen.io.vasp.inputs import Poscar
import numpy as np
import numpy.linalg as LA

class PoscarAdjuster:
    def __init__(self, poscar):
        self.poscar = poscar
    def modify_lc(self, lc):
        d = self.poscar.as_dict()
        A = d['structure']['lattice']['matrix']
        vacspc = A[-1,-1]
        lc0 = LA.norm(A[0])
        assert np.isclose(LA.norm(A[0]), LA.norm(A[1]))
        A = A / lc0 * lc
        A[-1,-1] = vacspc
        d['structure']['lattice']['matrix'] = A
        return (lc, Poscar.from_dict(d))

