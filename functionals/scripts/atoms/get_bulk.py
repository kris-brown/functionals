from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms # type: ignore

from ase.io import write,read # type: ignore

try:
    from bulk_enumerator.bulk import BULK # type: ignore
except ImportError:
    pass

from os      import remove,environ
from os.path import join
from random  import choices
from string  import ascii_lowercase
################################################################################
def get_bulk(atoms : 'Atoms'
            ,tol   : float  = 0.05
            ) -> 'BULK':
    """
    Create a BULK object from Ankit's library
    """
    # Constants
    #----------
    randroot = environ['DBGEN_STORAGE']
    suffix   = 'tmp_'+''.join(choices(ascii_lowercase,k=8))
    pth      = join(randroot,suffix)
    # Initialize
    #-----------
    assert tol <= 0.1 # requirement?
    b = BULK(tolerance=tol)
    poscar = ''
    # Get POSCAR String
    #------------------
    while poscar == '':
        write(pth,atoms,format='vasp')
        try:
            with open(pth,'r') as f:
                poscar = f.read()
        except IOError:
            pass # try repeatedly

    # Create bulk from POSCAR string
    #-------------------------------
    b.set_structure_from_file(poscar)

    if b.get_name()=='':
        raise ValueError('problem with getting bulk: \n'+poscar)
    try:
        remove(pth)
    except IOError:
        pass

    return b
