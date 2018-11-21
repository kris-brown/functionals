import numpy as np  # type: ignore
import ase          # type: ignore

def get_system_type(atoms : ase.Atoms) -> str:
    """
    Helpful docstring
    """
    ############################################################################
    # Parameters
    #------------
    cutoff    = 6  # A
    cell_abc  = list(map(np.linalg.norm,atoms.get_cell()[:3]))

    big_boy   = atoms.copy()*(2,2,2)
    xs,ys,zs  = zip(*big_boy.positions)
    sx,sy,sz  = [sorted(q) for q in [xs,ys,zs]]
    gap       = [max(np.diff(s)) for s in [sx,sy,sz]]
    pbc       = [gap[i] < cutoff for i in range(3)]

    if   all(pbc): return 'bulk'
    elif any(pbc): return 'surface'
    else:          return 'molecule'

if __name__=='__main__':
    import sys
    from dbgen.core.misc import anytraj
    a = anytraj(sys.argv[1])
    print(get_system_type(a))
