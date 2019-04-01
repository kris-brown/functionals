from typing import Callable as C, Tuple as T

def parse_incar(pth:str)->T[float,float,str,str,str,float,str,int,int,int,int,int,int,int,int,float,float,float,float,float,float,int]:

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
                magmom    = int)
    d = {x.strip().lower() : y.strip() for x,y in pairs}

    return tuple([maybe(f)(d.get(k)) for k,f in keys.items()]) # type: ignore
