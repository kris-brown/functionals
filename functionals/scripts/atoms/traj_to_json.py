import ase # type: ignore
import json
#from numpy.linalg import norm  # type: ignore

from dbgen.core.numeric import roundfloat
def traj_to_json(atoms : ase.Atoms) -> str:
    """
    Serialize an Atoms object in a human readable way
    """

    atomdata = []

    fixed_inds = []
    if atoms.constraints:
        for constraint in atoms.constraints:
            if isinstance(constraint,ase.constraints.FixAtoms):
                fixed_inds.extend(list(constraint.get_indices()))

    atoms.wrap()
    for a in atoms: atomdata.append({'number'       : int(a.number)
                                    ,'x'            : roundfloat(a.x)
                                    ,'y'            : roundfloat(a.y)
                                    ,'z'            : roundfloat(a.z)
                                    ,'magmom'       : roundfloat(a.magmom)
                                    ,'tag'          : int(a.tag)
                                    ,'constrained'  : int(a.index in fixed_inds)
                                    ,'index'        : int(a.index)})

    out = {'cell': [[roundfloat(x) for x in xx] for xx in atoms.get_cell().tolist()]
          ,'atomdata':atomdata}

    return json.dumps(out)
