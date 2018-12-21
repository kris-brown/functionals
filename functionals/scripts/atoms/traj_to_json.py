from ase import Atoms # type: ignore
from ase.constraints import FixAtoms # type: ignore

from json import dumps
#from numpy.linalg import norm  # type: ignore

from dbgen.utils.numeric import roundfloat
def traj_to_json(atoms : 'Atoms') -> str:
    """
    Serialize an Atoms object in a human readable way
    """
    def roundfloat(x: float) -> float:
        output = round(float(x), 3)
        return abs(output) if output == 0 else output

    atomdata = []

    fixed_inds = []
    if atoms.constraints:
        for constraint in atoms.constraints:
            if isinstance(constraint,FixAtoms):
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

    return dumps(out)
