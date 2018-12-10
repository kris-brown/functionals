# External modules
from ase.io     import read     # type: ignore
from ase        import Atoms    # type: ignore
from glob       import glob
from os.path    import getsize
# Internal modules
################################################################################
def anytraj(root : str) -> 'Atoms':
    """
    ASE IO read function - takes any traj it can find
    """
    trajs = glob('%s/*.traj'%root)
    for t in trajs:
        if getsize(t) > 10:
            atoms = read(t)
            return atoms
    raise ValueError('Get Traj could not find any traj in '+root)

if __name__=='__main__':
    import sys
    print(anytraj(sys.argv[1]))
