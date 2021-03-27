from __future__ import print_function
from typing import Optional, List, Dict
import glob
from ase import io, units
from ase.eos import EquationOfState
from sys import argv
import pylab
import os
import re
from prettytable import PrettyTable

'''
THIS IS FOR mCAML paper
'''


def gas_engs() -> Dict[str, float]:
    res = {}
    if os.path.exists('gas'):  # Kris
        for pth in os.listdir('gas'):
            if os.path.exists('gas/' + pth + '/OUTCAR'):
                res[pth] = get_eng('gas/' + pth + '/OUTCAR')
    elif os.path.exists('atoms'):  # Kai
        for pth in glob.glob('atoms/*_*/'):
            elem = pth.split('/')[1].split('_')[0]
            eng = get_eng(pth + '/OUTCAR')
            # Take the MAGMOM with lowest energy
            if res.get(elem, float('inf')) > eng:
                res[elem] = eng
    elif '/beefjobs/' in os.getcwd():
        gas = os.getcwd().replace('bulks', 'atoms')
        for pth in os.listdir(gas):
            res[pth] = get_eng(gas + '/' + pth + '/OUTCAR')
    else:
        raise ValueError()
    return res


def get_eng(pth: str) -> float:
    with open(pth, 'r') as f:
        outcar = f.read()
    pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
    match = re.findall(pat, outcar)

    if not match:
        raise IOError("Cannot parse energy %s - %s" % (match, pth))
    eng = float(match[-1])
    try:
        eng2 = io.read(pth).get_total_energy()
        assert eng2 == eng, (eng2, eng)
    except Exception:
        pass
    return float(eng)


def get_elems(name: str) -> List[str]:
    if len(name) == 1:
        return [name]
    elif str.isupper(name[1]):
        return [name[0], name[1:]]
    elif len(name) == 2:
        return [name]
    else:
        return [name[:2], name[2:]]


convsizefactors = {
    'fcc': 4.0,
    'bcc': 2.0,
    'diamond': 4.0,
    'rocksalt': 4.0,
    'zincblende': 4.0,
    'cesiumchloride': 1.0}
n_units = dict(fcc=1, bcc=1, diamond=2, rocksalt=1, zincblende=1,
               cesiumchloride=1)

scheffmats = dict(
    Ag=(2.93,112.60000000000001,4.063),
    Al=(3.36,81.0,4.019),
    AlAs=(3.7399999999999998,79.3,5.647),
    AlP=(4.209999999999999,87.0,5.448),
    AlSb=(3.26,59.5,6.121),
    Au=(3.79,182.9,4.061000000000001),
    Ba=(1.8900000000000001,9.3,5.003),
    BAs=(4.62,151.0,4.765),
    BN=(6.4399999999999995,388.5,3.593),
    BP=(4.9399999999999995,176.5,4.525),
    C=(7.1899999999999995,453.3,3.5530000000000004),
    Ca=(1.82,18.699999999999996,5.554),
    CdS=(2.7399999999999998,65.0,5.808),
    CdSe=(2.44,55.5,6.044),
    CdTe=(2.21,46.2,6.473000000000001),
    Cu=(3.4600000000000004,144.9,3.596),
    Fe=(4.24,177.9,2.858),
    GaAs=(3.2800000000000002,78.10000000000001,5.635),
    GaN=(4.4,204.0,4.5089999999999995),
    GaP=(3.5100000000000002,89.8,5.434),
    GaSb=(2.9699999999999998,56.8,6.076),
    Ge=(3.8300000000000005,77.9,5.646),
    HfC=(8.04,277.9,4.629),
    HfN=(7.9799999999999995,312.3,4.51),
    InAs=(3.02,58.4,6.029999999999999),
    InP=(3.4000000000000004,73.5,5.851999999999999),
    InSb=(2.77,47.199999999999996,6.461),
    Ir=(6.91,385.70000000000005,3.833),
    K=(0.93,3.7300000000000004,5.207999999999999),
    Li=(1.5899999999999999,13.4,3.4509999999999996),
    LiCl=(3.52,37.3,5.072),
    LiF=(4.34,75.4,3.973),
    MgO=(5.05,173.00000000000003,4.188999999999999),
    MgS=(3.96,81.0,5.188),
    Mo=(6.78,267.90000000000003,3.141),
    Na=(1.09,7.7,4.209),
    NaCl=(3.31,29.2,5.572),
    NaF=(3.89,54.099999999999994,4.582),
    Nb=(7.54,176.2,3.2929999999999997),
    NbC=(8.19,312.29999999999995,4.460999999999999),
    NbN=(7.47,300.5,4.3709999999999996),
    Ni=(4.4,190.7,3.5069999999999997),
    Pd=(3.87,197.4,3.876),
    Pt=(5.81,279.2,3.9130000000000003),
    Rb=(0.85,3.11,5.577),
    Rh=(5.72,271.6,3.794),
    Si=(4.56,100.3,5.420999999999999),
    SiC=(6.25,228.9,4.3469999999999995),
    Sn=(3.12,53.5,6.474),
    Sr=(1.71,12.5,6.04),
    Ta=(8.07,194.99999999999997,3.2990000000000004),
    TiC=(7.08,249.5,4.3180000000000005),
    TiN=(6.620000000000001,295.5,4.228),
    V=(5.279999999999999,158.79999999999998,3.023),
    VC=(6.88,307.7,4.149),
    VN=(6.2,281.29999999999995,4.122),
    W=(8.87,315.79999999999995,3.1599999999999997),
    ZnS=(3.13,78.10000000000001,5.392),
    ZnSe=(2.5999999999999996,65.9,5.6610000000000005),
    ZnTe=(2.39,53.099999999999994,6.09),
    ZrC=(7.87,229.9,4.686999999999999),
    ZrN=(7.46,220.3,4.574999999999999))

sol58mats = dict(
    AlN=(5.848, 206., 4.368),  # 5.845???
    LiH=(2.489, 40.1, 3.979),
    CaO=(5.520, None, 4.781),  # 5.534
    CoAl=(None, None, 2.854),
    FeAl=(None, None, 2.881),
    NiAl=(None, None, 2.881))
sol58lookup = dict(
    AlN='zincblende',
    LiH='rocksalt',
    CaO='rocksalt',
    CoAl='cesiumchloride',
    FeAl='cesiumchloride',
    NiAl='cesiumchloride')
mats = scheffmats # sol58mats

def mbfloat(x: str) -> Optional[float]:
    return None if not x else float(x)


kjmol_to_ev = units.kJ / units.mol


def main() -> None:
    gases = gas_engs()
    pt = PrettyTable(['mat', 'errCe', 'errBm', 'errLc'])
    ceErr, bmErr, lcErr = [], [], []
    tols = [2, 1, 3]
    s = sorted(glob.glob('*_*'))
    for i in range(len(s) - 1, -1, -1):
        if s[i].find('.') > -1:
            s.pop(i)

    for solid in s:
        elem, struc = solid.split('_')

        if struc not in n_units:
            elem = solid.replace('_', '')
            struc = sol58lookup.get(elem, '')

        if elem not in mats:
            print(solid, ' not in dataset')
            continue
        elems = get_elems(elem)
        volumes = []
        energies = []

        # Get gas energy
        gas_eng = sum(map(gases.__getitem__, elems))

        dirs = glob.glob(solid + '/dir*_0*')
        if not dirs:
            dirs = glob.glob(solid + '/eos/strain_*')
        if not dirs:  # KAI
            best = max([x.split('_')[-1] for x in glob.glob(solid + '/fiv*')])
            if elem == 'Na':
                best = str(int(best) - 1)
            dirs = glob.glob(solid + '/fivePointStencil_' + best + '/*')
        if not dirs:
            print('MISSING ', solid)
            continue
        for d in dirs:
            volumes.append(io.read(d + '/POSCAR').get_volume())
            energies.append(get_eng(d + '/OUTCAR'))

        eos = EquationOfState(volumes, energies)
        v0, E0, B = eos.fit()
        v0, E0 = sorted(zip(volumes, energies))[2]
        bulk_en = E0 / n_units[struc]  # bulkeng per chemical formula
        ce = (gas_eng - bulk_en) / len(elems)  # CE per atom
        B *= 1e24 / units.kJ  # GPa
        conventvol = v0 * convsizefactors[struc]
        latt = conventvol**(1. / 3.)
        exptCE, exptBM, exptLC = mats[elem]
        errCE = (ce - exptCE) if exptCE else ''
        errBM = B - exptBM if exptBM else ''
        errLC = latt - exptLC if exptLC else ''
        relBM = errBM / exptBM if errBM else ''
        relLC = errBM / exptBM if errBM else ''

        pt.add_row([elem] + [round(x, r) if x else '' for x, r in
                             zip([errCE, errBM, errLC], tols)])
        # print(errCE)
        if errCE:
            ceErr.append(float(errCE))
        if errBM:
            bmErr.append(float(errBM))
        if errLC:
            lcErr.append(float(errLC))
        if len(argv) > 1:
            eos.plot('%s/%s-eos.png' % (argv[1], elem))
            pylab.close()
    allErr = (ceErr, bmErr, lcErr)
    pt.add_row(['-'] * 4)
    pt.add_row(['mae'] + [str(round(sum(map(abs, x)) / len(x), r))
                          for x, r in zip(allErr, tols)])
    pt.add_row(['me'] + [str(round(sum(x) / len(x), r))
                         for x, r in zip(allErr, tols)])
    pt.add_row(['rmse'] + [str(round((sum(map(lambda y: y**2, x)
                                          ) / len(x))**0.5, r))
                           for x, r in zip(allErr, tols)])
    if mats == sol58mats:
        print(pt)
    else:
        print(len(ceErr))
        for e in ceErr: print(e)
        print('')
        for e in bmErr: print(e)
        print('')
        for e in lcErr: print(e)
        print('')


if __name__ == '__main__':
    main()
