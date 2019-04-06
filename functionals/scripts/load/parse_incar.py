from typing import Callable as C, Tuple as T

def parse_incar(pth:str)->dict:
    '''Turn path to INCAR into a dictionary'''
    with open(pth,'r') as f:
        pairs = [l.split('=') for l in f if len(l.split('=')) == 2]

    def boolean(x : str) -> bool:
        '''Parse a boolean from INCAR'''
        if   'true'  in x.lower(): return True
        elif 'false' in x.lower(): return False
        else:                      return bool(int(x))

    def maybe(f : C) -> C:
        return lambda x: x if x is None else f(x)

    keys = dict(encut     = float,
                sigma     = float,
                metagga   = str,
                gga       = str,
                prec      = str,
                ediff     = float,
                algo      = str,
                ismear    = int,
                npar      = int,
                nelm      = int,
                ispin     = int,
                ibrion    = int,
                lcharg    = boolean,
                lbeefens  = boolean,
                addgrid   = boolean,
                lasph     = boolean,
                lwave     = boolean,
                a11       = float,
                a12       = float,
                a13       = float,
                a14       = float,
                a15       = float,
                msb       = float,
                magmom    = float)
    d = {x.strip().lower() : y.strip() for x,y in pairs}

    return {k:maybe(f)(d.get(k)) for k,f in keys.items()} # type: ignore
