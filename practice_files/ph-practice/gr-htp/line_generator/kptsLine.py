from pymatgen.io.vasp.inputs import Kpoints # Generate KPOINTS line
from pymatgen.core.structure import Structure # Import structure obj from POSCAR
from pymatgen.symmetry.bandstructure import HighSymmKpath # Generate IBZ object

structure = Structure.from_file('POSCAR_unit')
print(structure)

ibz_fromStructure = HighSymmKpath(structure)

kpts_new = Kpoints.automatic_linemode(50, ibz_fromStructure)
print("New kpts:\n{}".format(kpts_new))