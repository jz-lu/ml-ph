from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, PH_FORCE_CONSTANTS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME,
    POSCAR_NAME
)
from ___constants_phonopy import POSCAR_UNIT_NAME
from ___constants_misc import ERR_PH_FORCE_SETS_NOT_MADE, ERR_PH_FORCE_CONSTS_NOT_MADE
from ___helpers_parsing import greet, succ, warn, err, update
from __directory_searchers import checkPath, findDirsinDir
from __dirModifications import build_dir
from pymatgen.io.vasp.inputs import Poscar
import os, phonopy, sys
import numpy as np; import numpy.linalg as LA
from hiphive import ForceConstants, ClusterSpace, ForceConstantPotential
from hiphive import enforce_rotational_sum_rules
from hiphive.utilities import extract_parameters
from ase.io import read
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

class PhonopyAPI:
    def __init__(self, ROOT, spname=SPOSCAR_NAME, pname=POSCAR_UNIT_NAME, ctype='twist'):
        greet("Initialized Phonopy API object. Loading phonopy objects...")
        warn("Assuming FORCE_SETS have already been extracted.")
        self.POSCAR_ERR_MSG = f"No unitcell POSCAR with name {pname} found"
        self.ROOT = checkPath(os.path.abspath(ROOT)); self.spname = spname; self.ctype = ctype
        self.nlayer = None; self.intra_list = None; self.nconfig = None; self.inter_list = None

        # Extract supercell dimensions, which may be different for 
        # pristine (intra) and configuration (inter) terms.
        d = checkPath(findDirsinDir(self.ROOT, MONOLAYER_DIR_NAME, searchType='start')[0])
        d = build_dir([d, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
        A0 = Poscar.from_file(d+pname).structure.lattice.matrix
        A1 = Poscar.from_file(d+spname).structure.lattice.matrix

        self.intra_sp_mat = np.round(A1 @ LA.inv(A0)) # diagonal matrix, diagonal is supercell dim
        d = checkPath(findDirsinDir(checkPath(self.ROOT + CONFIG_DIR_NAME), CONFIG_SUBDIR_NAME, searchType='start')[0])
        d = build_dir([d, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
        A2 = Poscar.from_file(d+spname).structure.lattice.matrix
        self.inter_sp_mat = np.round(A2 @ LA.inv(A0))
        print(f"Monolayer supercell size: {np.diag(self.intra_sp_mat)}")
        print(f"Configuration supercell size: {np.diag(self.inter_sp_mat)}")

        assert ctype in ['twist', 'config'], f"Unknown computation type {self.ctype}"
        if ctype == 'twist':
            self.nlayer, self.intra_list = self.__load_intra_ph_list(ROOT, pname)
            succ("Successfully loaded intralayer objects")
        self.nconfig, self.inter_list = self.__load_inter_ph_list(ROOT, pname)
        succ("Successfully loaded interlayer objects")

    # Load monolayer phonopy objects and apply rotational sum rule
    def __load_intra_ph_list(self, ROOT, pname=POSCAR_UNIT_NAME):
        ROOT = checkPath(ROOT); ph_list = []
        layers = sorted(findDirsinDir(ROOT, MONOLAYER_DIR_NAME, searchType='start'))
        layers.sort(key=len)
        assert len(layers) > 1, "Twist calculations require at least 2 layers"
        assert len(layers) <= 2, "Twist calculations for more than 2 layers not supported (yet)"
        cutoff = 4; dim = int(self.intra_sp_mat[0,0])
        for layer in layers:
            ph_dir = build_dir([ROOT, layer, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
            assert os.path.isdir(ph_dir), f"Directory {ph_dir} not found"
            os.chdir(ph_dir)
            assert os.path.isfile(pname), self.POSCAR_ERR_MSG
            assert os.path.isfile(PH_FORCE_CONSTANTS_NAME), ERR_PH_FORCE_CONSTS_NOT_MADE
            
            # supercell = read(SPOSCAR_NAME); prim = read(pname)
            # fcs_phonopy = ForceConstants.read_phonopy(supercell, PH_FORCE_CONSTANTS_NAME)
            # cs = ClusterSpace(prim, [cutoff])
            # parameters = extract_parameters(fcs_phonopy, cs)
            # enforced_parameters = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang'])
            # fcp_rot = ForceConstantPotential(cs, enforced_parameters)
            # fcs_hiphive_rot = fcp_rot.get_force_constants(supercell)
            # phonopy_prim = PhonopyAtoms(numbers=prim.numbers, 
            #                             positions=prim.positions, cell=prim.cell)
            # phon = Phonopy(phonopy_prim, 
            #                supercell_matrix=np.diag([dim, dim, 1]), 
            #                primitive_matrix=None)
            # phon.set_force_constants(fcs_hiphive_rot.get_fc_array(order=2))
            # ph_list.append(phon)
            
            ph_list.append(phonopy.load(unitcell_filename=pname, 
                                        supercell_matrix=np.diag([dim, dim, 1]), 
                                        primitive_matrix=np.eye(3),
                                        log_level=0)) # `log_level` is just the debug-paranoia level
            print(f"Loaded phonopy object from {ph_dir} and applied rotational sum rule")
        return len(layers), ph_list

    # Load configuration phonopy objects
    def __load_inter_ph_list(self, ROOT, pname=POSCAR_UNIT_NAME):
        if self.ctype == 'twist':
            ROOT = build_dir([ROOT, CONFIG_DIR_NAME])
        configs = sorted(findDirsinDir(ROOT, CONFIG_SUBDIR_NAME, searchType='start'))
        configs.sort(key=len)
        ph_list = []
        for config in configs:
            ph_dir = build_dir([ROOT, config, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
            assert os.path.isdir(ph_dir), f"Directory {ph_dir} not found"
            os.chdir(ph_dir)
            assert os.path.isfile(pname), self.POSCAR_ERR_MSG
            assert os.path.isfile(PH_FORCE_SETS_NAME), ERR_PH_FORCE_SETS_NOT_MADE
            ph_list.append(phonopy.load(unitcell_filename=pname, 
                                        supercell_matrix=np.diag(self.inter_sp_mat).astype(int), 
                                        primitive_matrix=np.eye(3),
                                        log_level=0)) # `log_level` is just the debug-paranoia level
        p = Poscar.from_file(build_dir([ROOT, configs[0]]) + POSCAR_NAME)
        self.bl_M = np.array(list(map(lambda x: x.atomic_mass, p.structure.species)))
        return len(configs), ph_list
    
    def nlayers(self):
        assert self.ctype == 'twist', f"Number of layers not defined for calculation type {self.ctype}"
        return self.nlayer
    def nconfigs(self):
        return self.nconfig
    def intra_ph_list(self):
        assert self.ctype == 'twist', f"Intralayer phonopy objects only built for twist calculations, not {self.ctype}"
        return self.intra_list
    def inter_ph_list(self):
        return self.inter_list
    def bl_masses(self):
        return self.bl_M


