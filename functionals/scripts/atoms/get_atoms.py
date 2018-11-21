from typing import List,Tuple

from json import loads

Ints = List[int]
Floats = List[float]
AtomTuple = Tuple[Ints,Ints,Floats,Floats,Floats,Ints,Floats,Ints]
################################################################################
def get_atoms(atomsjson : str)->AtomTuple:
    """
    Unpacks atom data from a json'd .traj file
    """
    ns,xs,ys,zs,ms,cs,ts,inds = [],[],[],[],[],[],[],[]
    atoms = loads(atomsjson)['atomdata']
    for a in atoms:
        ns.append(a['number']); xs.append(a['x']); ys.append(a['y'])
        zs.append(a['z']); ms.append(a['magmom']); inds.append(a['index']);
        cs.append(a['constrained']);ts.append(a['tag'])

    return ns,xs,ys,zs,cs,ms,ts,inds,ns # type: ignore
