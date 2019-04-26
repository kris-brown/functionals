from typing import Dict
from os.path import exists
import numpy  as np # type: ignore
from json   import dumps,loads
from dbgen.utils.parsing import parse_line


def parse_xc_vasp(stordir:str) -> str:
    """
    Finds the exchange correlation functional name, unless it is a BEEF functional
    in which case it returns a JSON of the coefficients
    """
    from os      import environ
    from os.path import join
    funroot = environ['FUNCTIONALS_ROOT']

    def readfile(pth:str)->str:
        with open(pth,'r') as f: return f.read()

    outcar = readfile(stordir+'/OUTCAR')

    parsed = parse_line(outcar,'GGA  ',0)
    if parsed is None:
        raise ValueError('malformed OUTCAR? no "GGA  " string found')
    else:
        gga = parsed.split()[2] # <name> <equals sign> <GGA name>
        if gga == 'PE':
            return 'PBE' # need check for if BEEF?
        elif gga == '--':
            beefcar = stordir+'/BEEFCAR'
            std     = join(funroot,'data/beef.json')
            arr_str = readfile(beefcar if exists(beefcar) else std)
            arr   = np.array(loads(arr_str))
            return dumps(arr.tolist())
        else:
            raise NotImplementedError('New functional in vasp? '+parsed)
