import json
import ase                              # type: ignore
import  pymatgen as pmg                 # type: ignore
import  pymatgen.symmetry.analyzer as pysym  # type: ignore

def get_pointgroup(atoms : ase.Atoms) -> str:
	"""
    uses pymatgen.symmetry.analyzer
	"""
	atoms.center()
	symbs = atoms.get_chemical_symbols()
	pos   = atoms.get_positions()
	mol   = pmg.Molecule(symbs,pos)
	return pysym.PointGroupAnalyzer(mol).sch_symbol
