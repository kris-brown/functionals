from typing import List as L, Tuple as T

from json import loads

################################################################################
def get_atoms(atomsjson : str)->T[L[int],L[float],L[float],L[float],L[bool],L[int],L[int],L[int]]:
    """
    Unpacks atom data from a json'd .traj file
    """
    ns,xs,ys,zs,ms,cs,ts,inds = [],[],[],[],[],[],[],[]
    atoms = loads(atomsjson)['atomdata']
    for a in atoms:
        ns.append(a['number']); xs.append(a['x']); ys.append(a['y'])
        zs.append(a['z']); ms.append(a['magmom']); inds.append(a['index']);
        cs.append(a['constrained']);ts.append(a['tag'])

    return ns,xs,ys,zs,cs,ms,ts,inds # type: ignore
