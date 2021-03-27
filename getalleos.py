from __future__ import print_function
from typing import Dict, List, Callable as C
import glob
import re
import os
from ase import io, units
from ase.eos import EquationOfState
from sys import argv
import pylab
'''
THIS IS FOR mBULK paper
'''
# ce_expt, lc_expt, bm_expt
scheffmats = dict(
    Ag=(2.970, 4.063000, 112.600),
    Al=(3.420, 4.019000, 81.000),
    AlAs=(3.820, 5.647000, 79.300),
    AlP=(4.310, 5.448000, 87.000),
    AlSb=(3.34, 6.121, 59.5),
    Au=(3.830, 4.061000, 182.9),
    Ba=(1.910, 5.003000, 9.300),
    BAs=(4.780, 4.765000, 151.000),
    BN=(6.760, 3.593000, 388.500),
    BP=(5.140, 4.525000, 176.500),
    C=(7.550, 3.553000, 453.300),
    Ca=(1.860, 5.554000, 18.700),
    CdS=(2.82, 5.808, 65),
    CdSe=(2.48, 6.044, 55.5),
    CdTe=(2.25, 6.473, 46.2),
    Cu=(3.520, 3.596000, 144.900),
    Fe=(4.320, 2.858000, 177.9),
    GaAs=(3.340, 5.635000, 78.100),
    GaN=(4.560, 4.509000, 204.000),
    GaP=(3.610, 5.434000, 89.800),
    GaSb=(3.03, 6.076, 56.8),
    Ge=(3.910, 5.646000, 77.900),
    HfC=(8.18, 4.629, 277.900),
    HfN=(8.08, 4.51, 312.3),
    InAs=(3.080, 6.030000, 58.400),
    InP=(3.460, 5.852000, 73.500),
    InSb=(2.81, 6.461, 47.2),
    Ir=(6.970, 3.833000, 385.7),
    K=(0.930, 5.208000, 3.73),
    Li=(1.670, 3.451000, 13.400),
    LiCl=(3.580, 5.072000, 37.300),
    LiF=(4.460, 3.973000, 75.400),
    MgO=(5.190, 4.189000, 173.000),
    MgS=(4.040, 5.188000, 81),
    Mo=(6.860, 3.141000, 267.9),
    Na=(1.130, 4.209000, 7.700),
    NaCl=(3.350, 5.572000, 29.200),
    NaF=(3.970, 4.582000, 54.100),
    Nb=(7.600, 3.293000, 176.2),
    NbC=(8.330, 4.461000, 312.3),
    NbN=(7.530, 4.371000, 300.5),
    Ni=(4.480, 3.507000, 190.7),
    Pd=(3.910, 3.876000, 197.400),
    Pt=(5.870, 3.913000, 279.2),
    Rb=(0.850, 5.577000, 3.11),
    Rh=(5.780, 3.794000, 271.600),
    Si=(4.700, 5.421000, 100.300),
    SiC=(6.470, 4.347000, 228.900),
    Sn=(3.160, 6.474000, 53.500),
    Sr=(1.730, 6.040000, 12.500),
    Ta=(8.110, 3.299000, 195),
    TiC=(7.240, 4.318000, 249.5),
    TiN=(6.760, 4.228000, 295.5),
    V=(5.340, 3.023000, 158.8),
    VC=(7.000, 4.149000, 307.7),
    VN=(6.300, 4.122000, 281.3),
    W=(8.930, 3.160000, 315.8),
    ZnS=(3.21, 5.392, 78.1),
    ZnSe=(2.66, 5.661, 65.9),
    ZnTe=(2.43, 6.09, 53.1),
    ZrC=(7.990, 4.687000, 229.9),
    ZrN=(7.580, 4.575000, 220.3)
)
othermats = dict(
    AgBr=(2.58, 5.78, None),
    AgCl=(2.78, 5.56, None),
    AgF=(2.95, 4.92, None),
    BaO=(5.1, 5.52, None),
    BaSe=(4.03, 6.6, None),
    CaSe=(4.01, 5.8, None),
    CdO=(3.21, 4.7, None),
    CaS=(4.81, 5.7, None),
    CoC=(5.69, 4.05, None),
    CoN=(4.53, 4.1, None),
    CrC=(5.8, 4.12, None),
    CrN=(5.14, 4.15, None),
    CsF=(3.66, 6.02, None),
    CsI=(2.75, 4.42, None),
    FeC=(5.67, 4.09, None),
    FeN=(4.59, 4.13, None),
    IrC=(6.84, 4.129, None),
    IrN=(5.13, 4.074, None),
    KBr=(3.08, 6.6, None),
    LaC=(3.08, 6.6, None),
    LaN=(6.27, 5.305, None),
    LiI=(2.78, 6, None),
    MnC=(5.14, 4.12, None),
    MnN=(4.08, 4.2, None),
    MnO=(4.75, 4.44, None),
    MnS=(4.01, 5.22, None),
    MoC=(7.22, 4.278, None),
    MoN=(6.2, 4.214, None),
    NiC=(5.65, 3.99, None),
    NiN=(4.48, 4.1, None),
    OsC=(7.36, 4.176, None),
    OsN=(5.61, 4.058, None),
    PdC=(4.03, 4.145, None),
    PdN=(5.36, 4.221, None),
    PtC=(6.34, 4.206, None),
    PtN=(4.63, 4.137, None),
    RbI=(2.71, 7, None),
    RhC=(6.23, 4.145, None),
    RhN=(4.78, 4.082, None),
    RuC=(6.73, 4.129, None),
    RuN=(5.02, 4.058, None),
    ScC=(6.37, 4.72, None),
    ScN=(6.72, 4.51, None),
    SeAs=(2.46, 5.48, None),
    TaC=(8.56, 4.457, None),
    TaN=(7.63, 4.34, None),
    WC=(8.25, 4.266, None),
    WN=(7.01, 4.202, None)
)
sol58mats = dict(
    AlN=(5.848, 4.368000, 206.000),  # 5.845???
    CaO=(5.520, 4.781000, None),  # 5.534
    CoAl=(None, 2.854000, None),
    FeAl=(None, 2.881000, None),
    LiH=(2.489, 3.979000, 40.100),
    NiAl=(None, 2.881000, None))

mats = scheffmats
print('USING SCHEFF MATS ONLY')
structs = dict(Ag='fcc',
               AgBr='rocksalt',
               AgCl='rocksalt',
               AgF='rocksalt',
               Al='fcc',
               AlAs='zincblende',
               AlN='zincblende',
               AlP='zincblende',
               AlSb='zincblende',
               Au='fcc',
               BAs='zincblende',
               BN='zincblende',
               BP='zincblende',
               Ba='bcc',
               BaO='rocksalt',
               BaSe='rocksalt',
               Be='hcp',
               C='diamond',
               Ca='fcc',
               CaO='rocksalt',
               CaS='rocksalt',
               CaSe='rocksalt',
               Cd='hcp',
               CdO='rocksalt',
               CdS='zincblende',
               CdSe='zincblende',
               CdTe='zincblende',
               Co='hcp',
               CoAl='cesiumchloride',
               CoC='rocksalt',
               CoN='rocksalt',
               CrC='rocksalt',
               CrN='rocksalt',
               CsF='rocksalt',
               CsI='cesiumchloride',
               Cu='fcc',
               Fe='bcc',
               FeAl='cesiumchloride',
               FeC='rocksalt',
               FeN='rocksalt',
               GaAs='zincblende',
               GaN='zincblende',
               GaP='zincblende',
               GaSb='zincblende',
               Ge='diamond',
               HfC='rocksalt',
               HfN='rocksalt',
               InAs='zincblende',
               InP='zincblende',
               InSb='zincblende',
               Ir='fcc',
               IrC='rocksalt',
               IrN='rocksalt',
               K='bcc',
               KBr='rocksalt',
               LaC='rocksalt',
               LaN='rocksalt',
               Li='bcc',
               LiCl='rocksalt',
               LiF='rocksalt',
               LiH='rocksalt',
               LiI='rocksalt',
               Mg='hcp',
               MgO='rocksalt',
               MgS='rocksalt',
               MnC='rocksalt',
               MnN='rocksalt',
               MnO='rocksalt',
               MnS='rocksalt',
               Mo='bcc',
               MoC='rocksalt',
               MoN='rocksalt',
               Na='bcc',
               NaCl='rocksalt',
               NaF='rocksalt',
               Nb='bcc',
               NbC='rocksalt',
               NbN='rocksalt',
               Ni='fcc',
               NiAl='cesiumchloride',
               NiC='rocksalt',
               NiN='rocksalt',
               Os='hcp',
               OsC='rocksalt',
               OsN='rocksalt',
               Pd='fcc',
               PdC='rocksalt',
               PdN='rocksalt',
               Pt='fcc',
               PtC='rocksalt',
               PtN='rocksalt',
               Rb='bcc',
               RbI='rocksalt',
               Rh='fcc',
               RhC='rocksalt',
               RhN='rocksalt',
               Ru='hcp',
               RuC='rocksalt',
               RuN='rocksalt',
               Sc='hcp',
               ScC='rocksalt',
               ScN='rocksalt',
               SeAs='rocksalt',
               Si='diamond',
               SiC='zincblende',
               Sn='diamond',
               Sr='fcc',
               Ta='bcc',
               TaC='rocksalt',
               TaN='rocksalt',
               Ti='hcp',
               TiC='rocksalt',
               TiN='rocksalt',
               V='bcc',
               VC='rocksalt',
               VN='rocksalt',
               W='bcc',
               WC='rocksalt',
               WN='rocksalt',
               Zn='hcp',
               ZnS='zincblende',
               ZnSe='zincblende',
               ZnTe='zincblende',
               Zr='hcp',
               ZrC='rocksalt',
               ZrN='rocksalt')
convsizefactors = {
    'fcc': 4.0,
    'bcc': 2.0,
    'diamond': 4.0,
    'rocksalt': 4.0,
    'zincblende': 4.0,
    'cesiumchloride': 1.0
}

n_units = dict(fcc=1, bcc=1, diamond=2, rocksalt=1, zincblende=1,
               cesiumchloride=1)


def stencil(engs: List[float], vols: List[float]) -> float:
    dx = vols[1] - vols[0]
    stenc = (-1 * engs[4] + 16 * engs[3] - 30 * engs[2] +
             16 * engs[1] - engs[0]) / (12 * dx**2)
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3

    bm = stenc * vols[2] * ev_a3_to_gpa  # APPROXIMATE
    return bm


def get_eng(pth: str) -> float:
    with open(pth, 'r') as f:
        outcar = f.read()
    pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
    match = re.findall(pat, outcar)

    if not match:
        raise IOError("Cannot parse energy %s - %s" % (match, pth))
    eng = float(match[-1])

    return float(eng)


def gas_eng(xc: str) -> C[[str], float]:
    gasdir = '/nfs/slac/g/suncatfs/ksb/beefjobs/atoms/%s/' % xc

    def f(elem: str) -> float:
        return get_eng(gasdir + elem + '/OUTCAR')
    return f


showCE, showBM, showLC = .3, .3, .3


def main(xc: str) -> None:
    make_figs, verbose = len(argv) > 2 and argv[2], len(argv) > 3

    bulkdir = '/nfs/slac/g/suncatfs/ksb/beefjobs/bulks/%s/' % xc

    errs: Dict[str, List[float]] = dict(ce=[], bm=[], bmstencil=[], lc=[])

    for name in mats.keys():
        if len(name) == 1:
            elems = [name]
        elif str.isupper(name[1]):
            elems = [name[0], name[1:]]
        elif len(name) == 2:
            elems = [name]
        else:
            elems = [name[:2], name[2:]]

        struc = structs[name]
        solid = '_'.join(elems)
        n_atom = len(elems)
        volumes, energies = [], []
        ce_expt, lc_expt, bm_expt = mats[name]
        dirs = sorted(glob.glob(bulkdir + '%s/eos/strain_*' % solid))
        if not dirs:
            print(bulkdir, solid, "NOT FOUND")
            continue
        for d in dirs:
            volumes.append(io.read(d + '/POSCAR').get_volume())
            try:
                energies.append(get_eng(d + '/OUTCAR'))
            except Exception as e:
                print(bulkdir, solid, e)
                continue
        if len(volumes) != len(energies):
            continue
        volumes, energies = zip(*sorted(zip(volumes, energies)))
        eos = EquationOfState(volumes, energies)
        v0, E0, B = eos.fit()
        B *= 1e24 / units.kJ  # GPa
        conventvol = v0 * convsizefactors[struc]
        latt = conventvol**(1. / 3.)

        if ce_expt:
            bulk_en = (E0 / n_units[struc])
            gas_en = sum(map(gas_eng(xc), elems))
            ce = (gas_en - bulk_en) / n_atom
            #
            errs['ce'].append(ce - ce_expt)
            if verbose or abs(errs['ce'][-1] / ce_expt) > showCE:
                print('CE ', solid, ce, ce_expt)
                # breakpoint()

        errs['lc'].append(latt - lc_expt)
        if verbose or abs(errs['lc'][-1] / lc_expt) > showLC:
            print('LC ', solid, latt, lc_expt)

        if bm_expt:
            errs['bmstencil'].append(abs(stencil(energies, volumes) - bm_expt))
            errs['bm'].append(B - bm_expt)
            if verbose or abs(errs['bm'][-1] / bm_expt) > showBM:
                print('BM ', solid, B, bm_expt)

        if make_figs:
            os.makedirs('eos_%s' % xc, exist_ok=True)
            eos.plot('eos_%s/%s-eos.png' % (xc, name))
            pylab.close()

    for met in ['ce', 'bm', 'bmstencil', 'lc']:
        if errs[met]:
            print(met, ' mae ', sum(map(
                lambda x: abs(x), errs[met])) / len(errs[met]))
        else:
            print('No ', met, 'data!')


if __name__ == '__main__':
    xc = argv[1]
    assert xc in ['pbe', 'ms2', 'pbesol', 'scan', 'mcaml', '5558', '7500',
                  '2751']
    main(xc)
