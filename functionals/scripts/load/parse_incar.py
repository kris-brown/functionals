from typing import Callable as C, Tuple as T


def parse_incar(pth: str) -> dict:
    '''Turn path to INCAR into a dictionary'''
    with open(pth, 'r') as f:
        pairs = [l.split('=') for l in f if len(l.split('=')) == 2]

    def boolean(x: str) -> bool:
        '''Parse a boolean from INCAR'''
        if 'true' in x.lower():
            return True
        elif 'false' in x.lower():
            return False
        else:
            return bool(int(x))

    def maybe(f: C) -> C:
        return lambda x: x if x is None else f(x)

    strs = ['metagga', 'gga', 'prec', 'magmom', 'ferwe', 'ferdo']
    floats = ['encut', 'ediff', 'sigma', 'nupdown']
    ints = ['ismear', 'npar', 'nelm', 'ispin', 'ibrion', 'icharg']
    bools = ['lcharg', 'lbeefens', 'addgrid', 'lasph', 'lwave', 'ldiag']
    keys = dict(**{k: str for k in strs}, **{k: float for k in floats},
                **{k: int for k in ints}, **{k: bool for k in bools})
    d = {x.strip().lower(): y.strip() for x, y in pairs}
    # make sure unique
    out = {k: maybe(f)(d.get(k)) for k, f in keys.items()}
    if out['magmom']:
        out['magmom'] = float(out['magmom'].replace('1*', '').split()[0])
    elif out['nupdown']:
        out['magmom'] = out['nupdown']
    return out
