# External modules
from ase.io     import read     # type: ignore
from glob       import glob
from os.path    import getsize
# Internal modules
from catalysis_model.scripts.Pure.Atoms.traj_to_json import traj_to_json
################################################################################
def anytraj(root : str) -> str:
    """
    ASE IO read function - takes any traj it can find
    """
    trajs = glob('%s/*.traj'%root)
    for t in trajs:
        if getsize(t) > 10:
            atoms = read(t)
            return traj_to_json(atoms)
    raise ValueError('Get Traj could not find any traj in '+root)

if __name__=='__main__':
    import sys
    print(anytraj(sys.argv[1]))
