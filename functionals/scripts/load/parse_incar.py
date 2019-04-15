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

    strs = ['metagga','gga','prec','magmom']
    floats = ['encut','ediff','sigma','a11','a12','a13','a14','a15','msb','nupdown']
    ints = ['ismear','npar','nelm','ispin','ibrion']
    bools = ['lcharg','lbeefens','addgrid','lasph','lwave']
    keys = dict(**{k:str for k in strs},**{k:float for k in floats},
                **{k:int for k in ints},**{k:bool  for k in bools})
    d = {x.strip().lower() : y.strip() for x,y in pairs}

    out = {k:maybe(f)(d.get(k)) for k,f in keys.items()} # type: ignore
    if   out['magmom']: out['magmom']=float(out['magmom'].split()[0])
    elif out['nupdown']: out['magmom'] = out['nupdown']
    return out
