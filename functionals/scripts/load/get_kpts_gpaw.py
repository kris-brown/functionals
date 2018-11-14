from typing import Dict,Tuple
import ast

from dbgen.core.parsing import parse_line

def get_kpts_gpaw(log:str) -> Tuple[int,int,int]:
    """
    docstring
    """
    parsed = parse_line(log,'k-points: ')
    if parsed is None:
        x,y,z=  (1,1,1) # default
    else:
        raw    = parsed.split(': ')[1].split()
        x,y,z = (int(raw[0]),int(raw[2]),int(raw[4]))

    assert all([x,y,z]), 'Kpoint parse error? raw=%s'%raw
    return x,y,z
