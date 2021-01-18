from __future__ import print_function
from typing import Dict, List
import glob
import re
import os
from ase import io, units
from ase.eos import EquationOfState
from sys import argv
import pylab


mats = dict(Ag_fcc=(2.970, 4.063000, 112.600),
            AlAs_zincblende=(3.820, 5.647000, 79.300),
            AlN_zincblende=(5.848, 4.368000, 206.000),
            AlP_zincblende=(4.310, 5.448000, 87.000),
            Al_fcc=(3.420, 4.019000, 81.000),
            Au_fcc=(3.830, 4.061000, None),
            BAs_zincblende=(4.780, 4.765000, 151.000),
            BN_zincblende=(6.760, 3.593000, 388.500),
            BP_zincblende=(5.140, 4.525000, 176.500),
            Ba_bcc=(1.910, 5.003000, 9.300),
            C_diamond=(7.550, 3.553000, 453.300),
            CaO_rocksalt=(5.520, 4.781000, None),
            Ca_fcc=(1.860, 5.554000, 18.700),
            CoAl_cesiumchloride=(None, 2.854000, None),
            Cu_fcc=(3.520, 3.596000, 144.900),
            FeAl_cesiumchloride=(None, 2.881000, None),
            Fe_bcc=(4.320, 2.858000, None),
            GaAs_zincblende=(3.340, 5.635000, 78.100),
            GaN_zincblende=(4.560, 4.509000, 204.000),
            GaP_zincblende=(3.610, 5.434000, 89.800),
            Ge_diamond=(3.910, 5.646000, 77.900),
            InAs_zincblende=(3.080, 6.030000, 58.400),
            InP_zincblende=(3.460, 5.852000, 73.500),
            Ir_fcc=(6.970, 3.833000, None),
            K_bcc=(0.930, 5.208000, None),
            LiCl_rocksalt=(3.580, 5.072000, 37.300),
            LiF_rocksalt=(4.460, 3.973000, 75.400),
            LiH_rocksalt=(2.489, 3.979000, 40.100),
            Li_bcc=(1.670, 3.451000, 13.400),
            MgO_rocksalt=(5.190, 4.189000, 173.000),
            MgS_rocksalt=(4.040, 5.188000, None),
            Mo_bcc=(6.860, 3.141000, None),
            NaCl_rocksalt=(3.350, 5.572000, 29.200),
            NaF_rocksalt=(3.970, 4.582000, 54.100),
            Na_bcc=(1.130, 4.209000, 7.700),
            NbC_rocksalt=(8.330, 4.461000, None),
            NbN_rocksalt=(7.530, 4.371000, None),
            Nb_bcc=(7.600, 3.293000, None),
            NiAl_cesiumchloride=(None, 2.881000, None),
            Ni_fcc=(4.480, 3.507000, None),
            Pd_fcc=(3.910, 3.876000, 197.400),
            Pt_fcc=(5.870, 3.913000, None),
            Rb_bcc=(0.850, 5.577000, None),
            Rh_fcc=(5.780, 3.794000, 271.600),
            SiC_zincblende=(6.470, 4.347000, 228.900),
            Si_diamond=(4.700, 5.421000, 100.300),
            Sn_diamond=(3.160, 6.474000, 53.500),
            Sr_fcc=(1.730, 6.040000, 12.500),
            Ta_bcc=(8.110, 3.299000, None),
            TiC_rocksalt=(7.240, 4.318000, None),
            TiN_rocksalt=(6.760, 4.228000, None),
            VC_rocksalt=(7.000, 4.149000, None),
            VN_rocksalt=(6.300, 4.122000, None),
            V_bcc=(5.340, 3.023000, None),
            W_bcc=(8.930, 3.160000, None),
            ZrC_rocksalt=(7.990, 4.687000, None),
            ZrN_rocksalt=(7.580, 4.575000, None))

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

tst = '/nfs/slac/g/suncatfs/ksb/tst/%s_%s'


def get_eng(pth: str) -> float:
    try:
        with open(pth, 'r') as f:
            outcar = f.read()
        pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
        match = re.findall(pat, outcar)

        if not match:
            raise IOError("Cannot parse energy %s - %s" % (match, pth))
        eng = float(match[-1])

        return float(eng)
    except Exception as e:
        raise e


def gas_eng(solid: str) -> float:
    elems = re.findall(r'[A-Z][a-z]?', solid)
    assert elems, solid
    tot = 0.

    for elem in elems:
        tot += get_eng('gas/' + elem + '/OUTCAR')
    return tot


def main(xc: str) -> None:

    s = sorted(glob.glob('*_*'))
    for i in range(len(s) - 1, -1, -1):
        if s[i].find('.') > -1:
            s.pop(i)

    errs: Dict[str, List[float]] = dict(ce=[], bm=[], lc=[])

    for solid in filter(lambda x: 'Pb' not in x, s):
        elem, struc = solid.split('_')
        n_atom = len(list(filter(str.isupper, elem)))
        volumes, energies = [], []
        ce_expt, lc_expt, bm_expt = mats[solid]

        xc_elem = tst % (xc, elem)
        dirr = xc_elem if os.path.exists(xc_elem) else solid
        dirs = sorted(glob.glob(dirr + '/dir*_0*'))

        for d in dirs:
            volumes.append(io.read(d + '/POSCAR').get_volume())
            energies.append(get_eng(d + '/OUTCAR'))
        eos = EquationOfState(volumes, energies)
        v0, E0, B = eos.fit()
        B *= 1e24 / units.kJ  # GPa
        conventvol = v0 * convsizefactors[struc]
        latt = conventvol**(1. / 3.)
        # print(elem, E0, latt, B)

        if ce_expt:
            ce = -((E0 / n_units[struc]) - gas_eng(elem)) / n_atom
            errs['ce'].append(ce - ce_expt)
            print(solid, ce, ce_expt)

        errs['lc'].append(latt - lc_expt)
        print(solid, latt, lc_expt)

        if bm_expt:
            errs['bm'].append(B - bm_expt)
            print(solid, B, bm_expt)

        if len(argv) > 2:
            eos.plot('%s/%s-eos.png' % (argv[2], elem))
            pylab.close()

    for met in ['ce', 'bm', 'lc']:
        print(met, ' mae ', sum(map(
            lambda x: abs(x), errs[met])) / len(errs[met]))


if __name__ == '__main__':
    xc = argv[1]
    assert xc in ['pbe', 'ms2', 'scan', 'mcaml']
    main(xc)
