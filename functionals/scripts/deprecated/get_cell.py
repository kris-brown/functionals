from typing import Tuple as T
from json import loads
###############################################################################
def get_cell(json_atoms : str)->T[float,float,float,float,float,float,float,float,float]:
    """
    Extracts information about a list of atoms (in json form) and prepares cell
    info for insertion into cell table
    """
    atoms = loads(json_atoms)
    [a,b,c] = atoms['cell']
    return  tuple(a+b+c) #type: ignore

if __name__=='__main__':
    import sys
    from catalysis_model.scripts.Pure.Atoms.traj_to_json import traj_to_json
    from ase.io import read  # type: ignore
    pth = sys.argv[1]
    print(get_cell(traj_to_json(read(pth))))
