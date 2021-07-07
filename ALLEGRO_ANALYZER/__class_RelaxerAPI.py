import numpy as np
import os
from ___constants_names import RELAX_CODE_PATH, RELAX_CODE_OUT
from __directory_searchers import checkPath
import numpy.linalg as LA

class RelaxerAPI:
    def __init__(self, theta, gridsz, outdir):
        assert isinstance(gridsz, int) and gridsz > 1, f"Invalid grid size {gridsz}: must be positive integer"
        assert theta > 0 and theta < 180, f"Invalid twist angle {theta}"
        assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        self.outdir = checkPath(os.path.abspath(outdir))
        self.outpath = self.outdir + RELAX_CODE_OUT
        self.gridsz = gridsz
        print("Starting relaxer program in Julia...")
        stream = os.popen(f"julia {RELAX_CODE_PATH} -d {theta} -N {gridsz} -o {self.outdir}")
        print(stream.read())
        assert os.path.isfile(self.outpath), f"Failed to find expected relaxer output file at {self.outpath}"
        print("Relaxer code finished.")
    def get_configs(self, cob):
        bprime_cart = np.load(self.outpath)
        assert bprime_cart.shape == (self.gridsz**2, 2), f"Invalid relaxation matrix shape (expected {(self.gridsz**2, 2)}):\n {bprime_cart}"
        bprime = np.round(((LA.inv(cob) @ bprime_cart.T).T + 1.000000001) % 1, 6)
        bprime = np.hstack((bprime, np.ones(self.gridsz**2).reshape(self.gridsz**2, 1)))
        return bprime

