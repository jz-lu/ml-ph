from ____exit_with_error import exit_with_error
from ___constants_config import *
from ___constants_names import POSCAR_NAME, POSCAR_CONFIG_NAMEPRE
from ___constants_misc import *
from ___constants_vasp import Z_LATTICE_SIZE, POSCAR_PRECISION, POSCAR_PREC_COMP_THRESHOLD, NO_RELAX_SELECTIVE_DYNAMICS_ARR, LAYER_RELAX_SELECTIVE_DYNAMICS_ARR
from __directory_searchers import findFilesInDir, checkPath
from pymatgen.io.vasp.inputs import Poscar # pylint: disable=import-error
import numpy as np # pylint: disable=import-error 
import copy

# Class Configuration
#    Stores the POSCAR and shift information for each displacement sampling.
#    Executes the strain-shift algorithm to prepare the composite layered material for relaxation.

class Configuration:
    # Import the POSCARs and check validity 
    def __init__(self, BASE_ROOT, poscars=None): # poscars: list of poscar objects imported with pmg
        self.BASE_ROOT = checkPath(BASE_ROOT)
        print('Configuration class instantiated. Initiating input files from %s'%(self.BASE_ROOT))

        # Import the given poscars if none are specified.
        if poscars == None:
            poscars = self.import_init_poscars()
        self.__poscars = tuple(poscars)
        print('POSCARs imported successfully.')


        # Get all the lattice information stored in self instance, and check validity of input.
        # self.__lattices: list of matrices of lattice basis vectors, each matrix for a POSCAR.
        print('Parsing lattice information...')
        self.__get_lattices()
        self.__get_normed_fixed_lattice()
        self.__get_lattice_constants()
        self.__check_lattice_consistency()
        self.__check_poscar_atoms()
        print('All basic consistency checks and verifications complete. Configuration object constructed.')

    # Import all the initial POSCARs and return them in a list of pmg poscar objects.
    # The input format for POSCARs in config sampling is 1 POSCAR per layer, minimum 2 layers.
    def import_init_poscars(self):
        # Every POSCAR in the root directory is imported.
        # NOTE: the fixed layer is the first one in alphabetic order
        poscar_names = findFilesInDir(self.BASE_ROOT, POSCAR_CONFIG_NAMEPRE, searchType='start')
        print("POSCAR files found:", poscar_names)
        print('Fixed layer set to POSCAR named ' + poscar_names[0])

        if len(poscar_names) < 2:
            exit_with_error(ERR_BAD_CONFIG_POSCAR)
        
        poscars = []
        for i in poscar_names:
            poscars.append(Poscar.from_file(self.BASE_ROOT + i))
        return poscars

    # Method to get the poscars safely.
    def get_init_poscars(self):
        return self.__poscars

    # For each imported POSCAR (corr. to each layer of the solid), get the lattice matrix, 
    # i.e. the 3 possibly un-normalized basis vectors concatenated into a column matrix.
    def __get_lattices(self):
        print('Extracting set of lattices from POSCARs...')
        # Get lattice matrices in an array.
        self.__lattices = []
        for i in self.__poscars:
            p_dict = i.as_dict()
            lattice = np.array(p_dict['structure']['lattice']['matrix'])
            self.__lattices.append(lattice)
        np.set_printoptions(precision=POSCAR_PRECISION)
        print('Lattices extracted. Ordered list:\n', self.__lattices)
        return self.__lattices

    # Extract the lattice constant of the POSCARs for straining.
    def __get_lattice_constants(self):
        print('Extracting set of lattice constants from lattices...')
        self.__lattice_constants = []
        for i in self.__lattices:
            diff = np.linalg.norm(i[0]) - np.linalg.norm(i[1])
            if diff > POSCAR_PREC_COMP_THRESHOLD: # first two lattice vectors need to have same norm always, float comp method
                exit_with_error(ERR_INCONSISTENT_NORMS)
            self.__lattice_constants.append(np.linalg.norm(i[0]))
            np.set_printoptions(precision=POSCAR_PRECISION)
        print('Lattice constants extracted. Ordered list:\n', self.__lattice_constants)
        return self.__lattice_constants

    # Get the normalized lattice matrix for the fixed (i.e. first) layer.
    def __get_normed_fixed_lattice(self):
        print('Retrieving normalized lattice basis of the fixed layer...')
        lattice = copy.deepcopy(self.__lattices[0])
        # Check that the in-plane lattice vectors have the same norm, otherwise everything is scaled wrong.
        diff = np.linalg.norm(lattice[0]) - np.linalg.norm(lattice[1]) 
        if diff != 0:
            print(WARN_LOW_INPLANE_PREC%(diff))
        if diff > POSCAR_PREC_COMP_THRESHOLD:
            print('Error in finding the normed fixed lattice matrix. In-plane lattice vectors have inconsistent norms.')
            exit_with_error(ERR_INCONSISTENT_NORMS)
        if lattice[2][2] != Z_LATTICE_SIZE: # Ensure that Z-plane lattice size is consistent with code's
            lattice[2][2] = Z_LATTICE_SIZE
            # exit_with_error(ERR_WRONG_ZVEC%(Z_LATTICE_SIZE))
        norm = np.linalg.norm(lattice[0])
        lattice[0] = lattice[0] / norm
        lattice[1] = lattice[1] / norm
        self.__normed_fixed_lattice = lattice
        print('Retrieval complete.')
        return lattice

    # Check if all the lattices are the same when they are normalized. If they are not, then
    # return error as strain-shift only works for layers of the same normalized lattice vectors.
    def __check_lattice_consistency(self):
        print('Verifying consistency of lattices (lattices must be identical up to constant multiples)...')
        # Loop over each lattice vector and check that they are the same as the fixed lattice vectors.
        for i in self.__lattices:
            for j in range(2):
                normed_vec = (i[j] / np.linalg.norm(i[j])) # Normalize
                abs_diff = np.linalg.norm(normed_vec - self.__normed_fixed_lattice[j])
                if abs_diff != 0:
                    print('Norm difference between fixed lattice and lattice vector %d of layer %d: %f'%(j+1, i+1, abs_diff))
                    print('Warning: precision of lattices in the nonfixed layers (wrt difference with fixed layer vectors) may not be sufficiently good for phonon calculations. Ensure all lattices and sublattices are of the same significant figures.')
                if abs_diff > POSCAR_PREC_COMP_THRESHOLD:
                    exit_with_error(ERR_INCONSISTENT_LATTICES)
            if i[2][2] != Z_LATTICE_SIZE: # Check that the Z-lattice initial spacing correct
                exit_with_error(ERR_WRONG_ZVEC%(Z_LATTICE_SIZE, i[2][2]))
        print('Lattice consistency verified.')
        return True

    # Check that the atoms in the same layer are coplanar.
    def __check_poscar_atoms(self):
        print('Checking coplanarity of atomic sublattice sites in POSCARs...')
        for i in range(len(self.__poscars)):
            mat = self.__poscars[i].structure.frac_coords
            for j in range(len(mat)):
                if mat[j][2] != mat[0][2]:
                    exit_with_error(ERR_ATOMS_OUT_OF_PLANE%(i+1, j+1))
        print('Atomic coplanarity validated.')

        # Return a list of selective dynamics bool arrays for interlayer relaxation
    
    # Return a SD indicator matrix constraining directions of relaxation.
    def __get_sd_matrix(self, num_fixed, num_nonfixed):
        # Interlayer relaxation allowed for all layers except the fixed layer.
        sd_mat = []
        for _ in range(num_fixed):
            sd_mat.append(NO_RELAX_SELECTIVE_DYNAMICS_ARR)
        for _ in range(num_nonfixed):
            sd_mat.append(LAYER_RELAX_SELECTIVE_DYNAMICS_ARR)
        return sd_mat

    # For each shift, construct a single POSCAR input file containing all the validated layers,
    # i.e. apply the strain-shift algorithm.
    def build_config_poscar(self, shift, init_interlayer_spacing):
        print('Building configuration POSCAR for sampling shift {}'.format(str(shift)))
        # Apply strain using the constant scalers array
        # Basically this is just scaling all the atoms in the other layers
        poscars = copy.deepcopy(self.__poscars) # So we don't delete everything in the class space

        # We need to strain (scale) everything 
        # so that it has the same lattice constant as the fixed layer.
        # lattice_constants = tuple(self.__lattice_constants)
        # lattice_scalers = []
        # for i in range(0, len(lattice_constants)):
        #     lattice_scalers.append(lattice_constants[0] / lattice_constants[i])
        # lattice_scalers = tuple(lattice_scalers) # no more modifying by accident

        bspace_structure = copy.deepcopy(poscars[0].structure)
        num_fixed_atoms = bspace_structure.num_sites # Number of atoms in fixed layer, needed for SD below
        num_nonfixed_atoms = 0 # Number of atoms we need to add to the first layer config space

        # Get a full structure object, except for SD, with strain-shift.
        for i in range(1, len(poscars)): # Don't modify the first (fixed) layer
            p = poscars[i]
            num_nonfixed_atoms += p.structure.num_sites
            n_at = p.structure.num_sites
            for _ in range(n_at): # Loop through each atom per layer sublattice
                at = p.structure.pop(0)
                # Shift it (no scaling is necessary since the lattice constant is already there!), 
                # modulo the unit cell torus which is just 1 in every coordinate in the lattice basis
                at.frac_coords = (at.frac_coords + shift) % 1
                
                # Since it is on a different layer we need to separate the layers in the z-coordinate
                at.frac_coords[2] = i * init_interlayer_spacing

                # Push it into the fxed lattice poscar object, which will be the b-space poscar
                bspace_structure.append(at.species, at.frac_coords)
            
        # SD: interlayer relaxation allowed for every 
        sd_mat = self.__get_sd_matrix(num_fixed_atoms, num_nonfixed_atoms)
        
        # Created a new fixed poscar with selective dynamics adjusted
        bspace_poscar = Poscar(bspace_structure, selective_dynamics=sd_mat)
        print('Build complete.')
        return bspace_poscar
        
    # Return set of poscar objects that each describe a particular shift in b-space
    # Return value: set of tuples (shift vector, poscar object)
    def build_config_poscar_set(self, shift_set, init_interlayer_spacing):
        print('Constructing the set of POSCARs now with shifts. Starting POSCAR builds...')
        configposcar_shift_tuple = []
        for i in shift_set:
            p = self.build_config_poscar(i, init_interlayer_spacing)
            configposcar_shift_tuple.append((i, p))

        self.config_space_poscar_set = configposcar_shift_tuple
        print('All shift poscar objects built.')
        return tuple(configposcar_shift_tuple)
        
    # Get the poscar of the fixed layer.
    def get_fixed_layer_poscar(self):
        return self.__poscars[0]
    
    @staticmethod
    # Returns a list of numpy row-vectors, each of which is a shift (expressed in arbitrary lattice basis).
    def sample_grid(grid=GRID_SAMPLE_LOW):
        sample_coord_sets = [] # All the sampling coordinates that we will zip together
        sample_points = [] # All possible combinations of the points in sample_coord_sets, with size grid[0]*grid[1]*grid[2]

        # Grid format validation.
        if len(grid) != 3 or grid[2] != 1: # Only shift in dimensions 1 annd 2 out of 3
            exit_with_error(ERR_INVALID_GRID)

        # Generate a set of shifts in each coordinate.
        for i in range(len(grid)):
            temp = []
            for sample_coord in range(grid[i]):
                temp.append(round((sample_coord / float(grid[i])) + 0.00000001, 6))
            sample_coord_sets.append(temp)
        
        # Construct a shift vector for every combination of the coordinate shifts.
        for i in sample_coord_sets[0]:
            for j in sample_coord_sets[1]:
                for k in sample_coord_sets[2]:
                    sample_point = (i, j, k)
                    sample_point = np.array(sample_point)
                    sample_points.append(sample_point)
        
        return tuple(sample_points)



# c = Configuration.sample_grid(GRID_SAMPLE_LOW)
# print(c)


# p = Poscar.from_file('POSCAR')
# q = Poscar.from_file('POSCAR')
# c = Configuration('.', [p,q])
# print(c.get_init_poscars())