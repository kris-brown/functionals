from typing import Tuple as T

######################################################
def get_kpts_vasp(stordir : str) -> T[int,int,int]:
    """
    docstring
    """
    with open(stordir+'/KPOINTS','r') as f:
        line = f.readlines()[3]
    x,y,z = tuple(map(int,line.split()))
    return x,y,z
