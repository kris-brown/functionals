import json
import numpy as np                      # type: ignore
import ase                              # type: ignore
from ase.constraints import FixAtoms    # type: ignore

def json_to_traj(raw_json:str)->ase.Atoms:
    """
    Inverse of traj_to_json
    """

    raw_atoms = json.loads(raw_json)

    atom_data = raw_atoms['atomdata']

    pos = np.array([[a[q] for a in atom_data] for q in ['x','y','z']]).T

    fix = FixAtoms([a['index'] for a in atom_data if a['constrained']])

    atoms = ase.Atoms(numbers     = [a['number'] for a in atom_data]
                     ,cell        = raw_atoms['cell']
                     ,positions   = pos
                     ,magmoms     = [a['magmom'] for a in atom_data]
                     ,tags        = [a['tag'] for a in atom_data]
                     ,constraint  = fix)
    return atoms
