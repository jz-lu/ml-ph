import numpy as np
import os
from ___constants_names import RELAX_CODE_PATH, RELAX_CODE_OUT
from __directory_searchers import checkPath

class RelaxerAPI:
    def __init__(self, theta, gridsz, outdir):
        assert isinstance(gridsz, int) and gridsz > 1, f"Invalid grid size {gridsz}"
        assert theta > 0 and theta < 180, f"Invalid twist angle {theta}"
        assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        self.outpath = checkPath(os.path.abspath(outdir)) + RELAX_CODE_OUT
        self.gridsz = gridsz
        print("Starting relaxer program in Julia...")
        stream = os.popen(f"julia {RELAX_CODE_PATH} -d {theta} -N {gridsz} -o {self.outpath}")
        print(stream.read())
        assert os.path.isfile(self.outpath), f"Failed to find expected relaxer output file at {self.outpath}"
        print("Relaxer code finished.")
    def get_u(self):
        u = np.load(self.outpath)
        assert u.shape == (self.gridsz**2, 2), f"Invalid relaxation matrix shape (expected {(self.gridsz**2, 2)}):\n {u}"
        return u

