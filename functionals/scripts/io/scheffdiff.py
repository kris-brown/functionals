from typing import Tuple as T, Optional as O, Dict as D


def scheffdiff(pth: str, name: str, calc: str, ce: O[float], bm: O[float],
               lc: O[float]) -> T[O[float], O[float], O[float]]:
    from csv import reader
    res_ce, res_bm, res_lc = (None, None, None)
    names = [
        'Ag', 'Al', 'AlAs', 'AlP', 'AlSb', 'Au', 'Ba', 'BAs', 'BN', 'BP',
        'C', 'Ca', 'CdS', 'CdSe', 'CdTe', 'Cu', 'Fe', 'GaAs', 'GaN', 'GaP',
        'GaSb', 'Ge', 'HfC', 'HfN', 'InAs', 'InP', 'InSb', 'Ir', 'K', 'Li',
        'LiCl', 'LiF', 'MgO', 'MgS', 'Mo', 'Na', 'NaCl', 'NaF', 'Nb', 'NbC',
        'NbN', 'Ni', 'Pd', 'Pt', 'Rb', 'Rh', 'Si', 'SiC', 'Sn', 'Sr', 'Ta',
        'TiC', 'TiN', 'V', 'VC', 'VN', 'W', 'ZnS', 'ZnSe', 'ZnTe', 'ZrC',
        'ZrN']
    res = dict()  # type: D[str,float]
    if calc in ['pbe', 'pbesol', 'scan'] and name in names:
        with open(pth, 'r') as f:
            r = reader(f)
            next(r)
            for (m, _, l1, l2, l3, xl, _, b1, b2, b3, xb, _,
                 c1, c2, c3, xc) in r:
                pbes = list(map(float, [c1, b1, l1]))
                pss = list(map(float, [c2, b2, l2]))
                scans = list(map(float, [c3, b3, l3]))
                xs = list(map(float, [xc, xb, xl]))
                for i, k in enumerate('cbl'):
                    res[m + k + 'pbe'] = (pbes[i] - xs[i])
                    res[m + k + 'pbesol'] = (pss[i] - xs[i])
                    res[m + k + 'scan'] = (scans[i] - xs[i])
        if ce is not None:
            sch = res[name + 'c' + calc]
            res_ce = float(ce) + sch
        if bm is not None:
            sch = res[name + 'b' + calc]
            res_bm = float(bm) - sch
        if lc is not None:
            sch = res[name + 'l' + calc]
            res_lc = float(lc) - sch
    return (res_ce, res_bm, res_lc)
