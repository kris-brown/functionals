from typing import Dict
from dbgen.utils.parsing import parse_line

def parse_pw_vasp(outcar:str) -> int:
    """
    Docstring
    """
    parsed  = parse_line(outcar,'ENCUT',0)
    if parsed is None:
        raise ValueError
    else:
        raw_num = parsed.split('=')[-1].split()[0]

    return int(round(float(raw_num)))
