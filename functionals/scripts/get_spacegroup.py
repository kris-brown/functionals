from typing import List,Tuple,TYPE_CHECKING
if TYPE_CHECKING:
	from ase import Atoms # type: ignore

from pymatgen.io.ase 			import AseAtomsAdaptor 	  # type: ignore
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer # type: ignore
from string import ascii_lowercase
def get_spacegroup(atoms : 'Atoms')->int: #Tuple[int,List[int],List[int]]:
	"""
    uses pymatgen SpacegroupAnalyzer
	"""
	pmg = AseAtomsAdaptor().get_structure(atoms)
	sd  = SpacegroupAnalyzer(pmg,symprec=0.1).get_symmetry_dataset()
	sg   = sd['number']
	#wy   = [ascii_lowercase.index(x) for x in sd['wyckoffs']]

	return sg # (sg,list(range(len(atoms))),wy)
