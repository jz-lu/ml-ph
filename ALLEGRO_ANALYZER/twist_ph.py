# Compute the phonon modes for a twisted cell
import numpy as np
from phonopy import Phonopy
import phonopy
from dm import TwistedDM, InterlayerDM, MonolayerDM
from bzsampler import BZSampler
from __directory_searchers import checkPath
from ___constants_names import SPOSCAR_NAME

# Load phonopy objects from all of the relevant directories
def load_ph_list(ROOT, spname=SPOSCAR_NAME):
    # From the monolayer directories, load the SPOSCAR for each layer
    ROOT = checkPath(ROOT) # TODO loop over a subdirectory for layer names as well (do the calculation part first!)
    monolayer_list = [phonopy.load(supercell_filename=ROOT+spname, force_constants_filename=)]
    return

# Sample moire G vectors (GM) and k-points along IBZ

