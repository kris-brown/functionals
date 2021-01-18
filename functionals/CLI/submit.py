from typing import Set, Dict, Any, Tuple, List, Optional
import argparse
import collections
import copy
import csv
import json
import math
import os
import os.path
import subprocess
import ase
import ase.data
import ase.io
import ase.eos
import numpy as np
import functionals.templates as temp
import scipy.optimize as opt
from ase.calculators.calculator import kptdensity2monkhorstpack

########
# DATA #
########


def str2bool(v: str) -> bool:
    return v.lower() in ("yes", "true", "t", "1")


gt = temp.jinja_env.get_template
ppth = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/potpaw_PBE'
root = '/nfs/slac/g/suncatfs/ksb/beefjobs/'
exe_root = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/'\
           'olsherlock/vasp.5.4.1beefcar/bin/vasp_'

data_root = '/' + os.path.join(*__file__.split('/')[:-3]) + '/data'
struct_root = data_root + '/structures'
vol_csv = data_root + '/opt_vol.csv'
beef_root = data_root + '/beefs'

beefs = [x[:-5] for x in os.listdir(beef_root)]

atommag, n_electrons = {}, {}
Mat = collections.namedtuple(
    'Mat', ['name', 'struct', 'lat', 'bm', 'mag', 'ca_rat'])

with open(data_root + '/elements.csv', 'r') as f:
    r = csv.reader(f)
    next(r)
    for e, n, m in r:
        atommag[e], n_electrons[e] = int(m), int(n)


with open(data_root + '/expt.csv', 'r') as f:
    r = csv.reader(f)
    has_data = {x: (bool(a), bool(b), bool(c), bool(d))
                for x, a, b, c, d in r}

with open(data_root + '/initialguess.csv', 'r') as f:
    r = csv.reader(f)
    next(r)
    matdata = {k: Mat(k, s, float(l), float(b), float(m or '0'),
                      float(ca or '0'))
               for k, s, l, b, m, ca in r
               if k in has_data}

mGGAincar = dict(NELMETA1=10, NELMETA2=10)
xcdict = dict(pbe=dict(gga='PE'), scan=dict(metagga='SCAN', **mGGAincar),
              pbesol=dict(gga='PS'),
              )  # type: Dict[str,Dict[str,Any]]
xcdict.update(**{beef: dict(metagga='BF', lbeefens=True, **mGGAincar)
                 for beef in beefs if beef not in xcdict})


def only_ce(m: Mat) -> bool:
    ce, bm, lat, _ = has_data[m.name]
    return ce and not (bm or lat)


def set_vol(a: ase.Atoms, v: float) -> None:
    r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
        a.get_cell() * x).volume - v, x0=0.5, x1=2)
    a.set_cell(
        a.get_cell() * r.root, scale_atoms=True)


##############################################################################
# HELPER FUNCTIONS #
n_dict = dict(fcc=4, bcc=2, hcp=2, cesiumchloride=2, diamond=8,
              zincblende=8, rocksalt=8, **{"rocksalt-prim": 2})


def read_outcar(pth: str) -> ase.Atoms:
    '''Zr gets mangled in CONTCAR which screws up OUTCAR parsing'''
    os.system("sed -i 's/ r/Zr/g' " + pth + '/CONTCAR')
    return ase.io.read(pth + '/OUTCAR')


def mk_traj(m: Mat) -> ase.Atoms:
    '''Create an Atoms object from name, structure, and lattice constant.'''
    c = [[1., 0, 0], [0, 1, 0], [0, 0, 1]]

    if m.struct == 'bcc':
        n = 1
        p = [[0., 0, 0]]
        c = [[-1 / 2, 1 / 2, 1 / 2], [1 / 2, -1 / 2, 1 / 2],
             [1 / 2, 1 / 2, -1 / 2]]
    elif m.struct == 'cesiumchloride':
        n = 1
        p = [[0, 0, 0], [1 / 2, 1 / 2, 1 / 2]]
    elif m.struct in ['hcp', 'wurtzite']:
        n = 2
        p = ([[0, 0, 0], [1 / 3, 2 / 3, 1 / 2]] if m.struct == 'hcp' else
             [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.376],
              [2 / 3, 1 / 3, 1 / 2], [2 / 3, 1 / 3, 0.876]])
        c = [[1, 0, 0], [-0.5, 3**0.5 / 2, 0], [0, 0, m.ca_rat]]  # DIFFERENT
    elif m.struct in ['zincblende', 'diamond']:
        n = 1 if m.struct == 'zincblende' else 2
        p = [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]
        c = [[0., 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
    elif m.struct in ['rocksalt', 'fcc']:
        n = 1
        p = [[0, 0, 0], [1 / 2, 1 / 2, 1 / 2]]
        c = [[0., 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
        if m.struct == 'fcc':
            p = p[:1]
    else:
        raise NotImplementedError(m.struct)
    return ase.Atoms(m.name * n, cell=np.array(c) * m.lat, scaled_positions=p,
                     pbc=[1, 1, 1],)


def binary(m: Mat) -> bool:
    return len(list(filter(str.isupper, m.name))) == 2


def snake(name: str) -> str:
    """Snake case materials with multiple elements"""
    newmat = name[0]
    for char in name[1:]:
        if char.isupper():
            newmat += "_" + char
        else:
            newmat += char
    return newmat


def onezero(one: int) -> str:
    return ('%d*1.0' % one if one else '') + ' %d*0.0' % (16 - one)  # FERWE/DO


def get_strain(mat: Mat, strain_int: int) -> float:
    '''Gives the *lattice* strain, in order to get a *volume* strain.'''
    spacing = float(0.25 / (mat.bm)**(1. / 3.))
    if mat in ['K', 'Ca', 'Rb', 'Sr']:
        spacing * 3
    return (1 + strain_int * spacing)**(1 / 3)


def current_jobs() -> Set[str]:
    cmd = 'bjobs -o "job_name" -noheader'
    return set(subprocess.check_output(cmd, encoding='UTF-8',
                                       shell=True).split())


def readfile(pth: str) -> str:
    """Convert filepath to string of contents."""
    with open(pth, 'r') as f:
        return f.read()


def done(spth: str) -> str:
    '''Identify if a job has already successfully completed.'''
    ozpth, outpth, rpth = [os.path.join(spth, f)
                           for f in ['OSZICAR', 'OUTCAR', 'run.log']]

    if not os.path.exists(ozpth):
        return 'NO OSZICAR'
    elif not os.path.exists(outpth):
        return 'NO OUTCAR'
    elif len(readfile(ozpth).split('\n')) > 800:
        return 'OSZICAR TOO LONG - UNCONVERGED?'
    elif 'General timing' not in readfile(outpth):
        return 'INCOMPLETE OUTCAR'
    elif 'fatal error' in readfile(rpth):
        return 'ERROR in run.log'
    return ''


def mk_incar(d: Dict[str, Any], pth: str) -> None:
    """Convert dictionary to INCAR file."""
    def fmt(val: Any) -> str:
        if isinstance(val, bool):
            return '.TRUE.' if val else '.FALSE.'
        else:
            return str(val)

    with open(pth, 'w') as f:
        lines = ['{} = {}'.format(k.upper(), fmt(v)) for k, v in d.items()]
        f.write('\n'.join(lines) + '\n')


def potpos(pth: str, atoms: ase.Atoms) -> None:
    poscar, potcar = [os.path.join(pth, x) for x in ['POSCAR', 'POTCAR']]
    ase.io.write(poscar, atoms)
    elems = subprocess.check_output('head -n 1 ' + poscar, encoding='UTF-8',
                                    shell=True).split()
    potcmd = ('cat ' + ' '.join(['%s/%s/POTCAR' % (ppth, x) for x in elems]
                                ) + ' > ' + potcar)
    os.system(potcmd)


def attempt(retry: bool, pth: str, curr: Set[str]) -> bool:
    notdone = bool(done(pth))
    notcurr = os.path.join(pth, 'subVASP.sh') not in curr
    return retry or (notdone and notcurr)


def magdic(mat: Mat) -> Dict[str, str]:
    if mat.mag:
        n = n_dict[mat.struct]
        b = binary(mat)
        mm = ' '.join(['0' if (b and (i % 2)) else str(mat.mag)
                       for i in range(n)])
        return {'ispin': "2", 'magmom': mm}
    else:
        return dict(ispin='1')


def mk_kpts(a: ase.Atoms) -> Tuple[int, int, int]:
    x, y, z = kptdensity2monkhorstpack(a, 6)
    return int(x), int(y), int(z)
##############################################################################
# MAIN FUNCTIONS #


def submit_atom(name: str, time: int, xc: str, retry: bool,
                incar: Dict[str, Any],
                curr: Set[str], orbitals: Optional[bool], sunc: str,
                beef: Optional[str]) -> None:
    assert sunc == '3'
    elem = ase.data.chemical_symbols.index(name)
    atoms = ase.Atoms(numbers=[elem], positions=[[1., 2., 3.]], pbc=[0, 0, 0],
                      cell=[[13., 2., 1.], [2., 14, 1.], [5., 1, 15.]])
    magmom, elec = atommag[name], n_electrons[name]
    nb = max(16, int(math.ceil(0.6 * elec + magmom
                               if magmom else elec / 2 + 0.5 * 1)))
    # ferwedo = orbitals if (orbitals is not None) else (elem > 56)
    # if ferwedo and magmom and name not in ['Pt', 'Au']:
    #     nup, ndn = (elec + magmom) // 2, (elec - magmom) // 2
    #     magdic = dict(ldiag=False, ismear=-2, ferwe=onezero(nup),
    #                   ispin=2, ferdo=onezero(ndn))
    # else:
    magdic = dict(sigma=0.0001, nupdown=magmom, ispin=2 if magmom else 1)

    magcar = dict(
        isym=-1, amix=0.1 if xc == 'scan' else 0.2,
        bmix=0.0001, nbands=nb, **magdic)

    # Submit
    pth = os.path.join(root, 'atoms', xc, name)
    os.makedirs(pth, exist_ok=True)
    if attempt(retry, pth, curr):
        setup(pth=pth, atoms=atoms, incar={**incar, **magcar,
                                           **dict(sigma=1e-4)},
              kpts=(1, 1, 1), beef=beef)
        bash = temp.jinja_env.get_template(
            'subVASP.jinja').render(gam='_gam', rm=True)
        _sub(pth, bash, time, sunc)


def phonon(mat: Mat, time: int, xc: str,
           retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
           beef: Optional[str]) -> None:
    raise NotImplementedError


def latopt(mat: Mat, time: int, xc: str,
           retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
           beef: Optional[str]) -> None:
    new = dict(ediff=1e-5, nsw=500, ediffg=-0.001, addgrid=True, algo='N',
               ismear=0, ivdw=0, lreal=False, ibrion=2, sigma=0.01,
               isif=3 if mat.struct in ['hcp', 'wurtzite'] else 7)
    atoms = mk_traj(mat)
    pth = os.path.join(root, 'bulks', xc, snake(mat.name), 'latopt')
    os.makedirs(pth, exist_ok=True)
    if attempt(retry, pth, curr):
        setup(pth=pth, atoms=atoms, incar={**incar, **new, **magdic(mat)},
              kpts=mk_kpts(atoms), beef=beef)
        bash = temp.jinja_env.get_template('subVASP.jinja'
                                           ).render(gam='', rm=False)
        _sub(pth, bash, time, sunc)


def eos(mat: Mat, time: int, xc: str,
        retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
        beef: Optional[str], second: bool) -> None:
    if second and only_ce(mat):
        return print('No EOS round 2 on %s without bm/lat data' % mat.name)
    matpth = os.path.join(root, 'bulks', xc, snake(mat.name),
                          'eos' + ('2' if second else ''))
    os.makedirs(matpth, exist_ok=True)

    if not second:
        outpth = matpth[:-3] + 'latopt'
        if done(outpth):
            return print('Cannot do EOS for %s - ' % matpth + done(outpth))
        try:
            orig_atoms = read_outcar(outpth)
            strain_vols = [ase.cell.Cell(orig_atoms.get_cell() *
                           get_strain(mat, i)).volume for i in range(-2, 3)]

        except Exception:
            return print('Cannot EOS for %s - ase parse error OUTCAR' % matpth)
    else:
        orig_atoms = copy.deepcopy(mk_traj(mat))
        output = matpth[:-1] + '/strain_%d'
        vols, engs = [], []
        for strain_i in [-2, -1, 0, 1, 2]:
            atoms = read_outcar((output % strain_i))
            vols.append(atoms.get_volume())
            engs.append(atoms.get_potential_energy())
        dx = vols[1] - vols[0]

        for i, j in [(0, 1), (1, 2), (3, 2), (4, 3)]:
            if engs[i] < engs[j]:
                print('BAD ENG ORDER: ',
                      ','.join(['%.2f' % x for x in engs]))
                return print('VOL ORDER: ',
                             ','.join(['%.2f' % x for x in vols]))
            if abs(abs(vols[i] - vols[j]) - dx) > .01:
                return print("UNEQUAL SPACING",
                             ','.join(['%.2f' % x for x in vols]))

        eos = ase.eos.EquationOfState(vols, engs)
        optvol, _, eosbm = eos.fit()  # type: Tuple[float,float,float]

        stencil = (-1 * engs[4] + 16 * engs[3] - 30 * engs[2] +
                   16 * engs[1] - engs[0]) / (12 * dx**2)
        ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3

        bulkmod = stencil * vols[2] * ev_a3_to_gpa  # APPROXIMATE
        eosbm = eosbm / ase.units.kJ * 1.0e24
        err = abs(bulkmod - eosbm) / eosbm
        # if err > 0.20:
        #     return print(xc, mat.name, " Irregular!", bulkmod, eosbm)

        r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
            orig_atoms.get_cell() * x).volume - optvol, x0=0.5, x1=2)
        orig_atoms.set_cell(
            orig_atoms.get_cell() * r.root, scale_atoms=True)
        v = orig_atoms.get_volume()
        assert round(v, 5) == round(optvol, 5), (output, v, optvol)

        shift = vols[2] - optvol
        strain_vols = [v - shift for v in vols]

    # Bulk-specific INCAR settings
    magcar = dict(sigma=0.01, ismear=0, **magdic(mat))

    dirs = []

    for svol, strain_i in zip(strain_vols, range(-2, 3)):
        # Strained Traj
        atoms = orig_atoms.copy()
        set_vol(atoms, svol)
        # Determine whether to submit again or not
        pth = os.path.join(matpth, 'strain_%d' % strain_i)
        os.makedirs(pth, exist_ok=True)
        if done(pth) != '':
            setup(pth=pth, atoms=atoms, incar={**incar, **magcar},
                  kpts=mk_kpts(atoms), beef=beef)
            dirs.append(pth)

    if dirs and attempt(retry, matpth, curr):

        bash = gt('subAll.jinja').render(gam='', dirs=dirs)
        _sub(matpth, bash, time, sunc)


def _sub(pth: str, bash: str, time: int, sunc: str) -> None:
    bashpth = os.path.join(pth, 'subVASP.sh')
    with open(bashpth, 'w') as g:
        g.write(bash)
    os.chdir(pth)
    os.system('chmod 755 ' + bashpth)
    args = [16 if sunc else 8, time, sunc, bashpth]
    cmd = 'bsub -n {} -W{}:09 -q suncat{} {}'.format(*args)
    os.system(cmd)


def setup(pth: str, atoms: ase.Atoms, incar: Dict[str, Any],
          kpts: Tuple[int, int, int], beef: Optional[str]) -> None:
    '''General setup of directory for atomic or bulk (at a strain).'''

    # Write INCAR
    mk_incar(incar, os.path.join(pth, 'INCAR'))

    # Write BEEFCAR
    if beef:
        with open(os.path.join(pth, 'BEEFCAR'), 'w') as f:
            f.write(beef)

    # KPOINTS
    with open(os.path.join(pth, 'KPOINTS'), 'w') as g:
        g.write(gt('KPOINTS.jinja').render(kpts=kpts))

    # Write POSCAR and POTCAR
    potpos(pth, atoms)


def restrain(pth: str, xc: str, strains: List[float]) -> None:
    dirs = []
    for i, strain in zip(range(-2, 3), strains):
        spth = os.path.join(pth, 'strain_%d' % i)
        dirs.append(spth)
        atoms = ase.io.read(spth + '/POSCAR')

        res = opt.root_scalar(
            f=lambda x: ase.atoms.Cell(atoms.get_cell() * x).volume - strain,
            x0=0.5, x1=1.5)

        atoms.set_cell(atoms.get_cell() * res.root, scale_atoms=True)
        ase.io.write(spth + '/POSCAR', atoms)
        for x in ['OUTCAR', 'WAVECAR', 'OSZICAR', 'EIGENVAL', 'REPORT']:
            try:
                os.remove(spth + '/' + x)
            except OSError:
                pass
    bash = gt('subAll.jinja').render(gam='', dirs=dirs)
    _sub(pth, bash, 10, '3')


#############################################################################
parser = argparse.ArgumentParser(description='Submit some jobs',
                                 allow_abbrev=True)
parser.add_argument('--time', default=10, type=int, help='Walltime')
parser.add_argument('--orb', default=None, type=str2bool,
                    help='Use explicit orb occupations')
parser.add_argument('--retry', default=False, type=bool, help='Redo jobs')
parser.add_argument('--sunc', default='', type=str, help='Which cluster',
                    choices=['', '3'])
parser.add_argument('--type', help='What type of bulk operation',
                    choices=['latopt', 'eos', 'phonon'])

parser.add_argument('--mat', default='', help='Material name',
                    type=lambda x: (list(matdata.values())
                                    if x == 'all' else
                                    [matdata[z] for z in x.split()]))
parser.add_argument('--elems', default='', help='"all" or list of pos ints',
                    type=lambda x: (list(atommag.keys())
                                    if x == 'all' else x.split()))
parser.add_argument('--xc', help='Copies xc into the working dir as BEEFCAR',
                    type=lambda x: (['pbe', 'pbesol', 'scan', 'ms2']
                                    if x == 'all' else x.split()))

parser.add_argument('--strains', default='', type=str,
                    help='Redo particular eos1 with modified volumes')

typedict = dict(latopt=latopt, eos=eos, phonon=phonon)  # type: Any


def main() -> None:
    '''Submit either bulk or element singlepoint calculations.'''
    # Get/validate command line args
    args = parser.parse_args()
    for xc in args.xc:
        if xc not in ['pbe', 'scan', 'pbesol']:
            assert xc in beefs
            with open(os.path.join(beef_root, xc + '.json'), 'r') as f:
                beef = ' '.join(map(str, json.load(f) + [0., 1.]))
        else:
            beef = None  # type: ignore
        assert args.time > 0
        curr = current_jobs()

        incar = dict(encut=900., ediff=1e-6, algo='N', istart=2,
                     addgrid=True, ibrion=-1, npar=1 if args.mat else 4,
                     lasph=True, nelm=800, lorbit=10,
                     lcharg=False, lwave=True, prec='ACCURATE', **xcdict[xc])

        common = dict(incar=incar, time=args.time, xc=xc, curr=curr,
                      retry=args.retry, sunc=args.sunc, beef=beef)
        if args.mat:
            assert (not args.elems) and args.type, args
            if args.strains:
                [m], [xc] = args.mat, args.xc
                eos = os.path.join(root, 'bulks', xc, snake(m.name), 'eos')
                eos2 = eos + '2'
                assert not os.path.exists(eos2) or not os.listdir(eos2), eos2
                assert args.strains.count(' ') == 4
                restrain(eos, xc, list(map(int, args.strains.split())))
                return None
            for m in args.mat:
                if args.type == 'eos':
                    matpth = os.path.join(root, 'bulks', xc, snake(m.name),
                                          'eos/')
                    second = os.path.exists(matpth)
                    if second:
                        for strain_i in [-2, -1, 0, 1, 2]:
                            strainpth = matpth + ('strain_%d') % strain_i
                            if done(strainpth):
                                # print('Cant do eos2 b/c incomplete %s - %s' %
                                #       (strainpth, done(strainpth)))
                                second = False
                                break
                    sdic = dict(second=second)
                else:
                    sdic = dict()
                typedict[args.type](mat=m, **sdic, **common)
        else:
            assert args.elems
            for e in args.elems:
                submit_atom(name=e, orbitals=args.orb, **common)


if __name__ == '__main__':
    main()

# HCP: Be Cd Co Hf In Mg Os Re Ru Sc Tc Ti Tl Y Zn Zr
# Ba_O Be Ca Ca_O Ca_S Ca_Se Cd Co Cr_C Cr_N Cs_F Fe_N Ga_N In_Sb Ir_C Ir_N K K_Br La_C La_N Li Mg Mg_O Mg_S Mn_C Mn_N Mn_O Mo Na Na_F Nb Nb_N Ni Ni_C Os_C Os_N Pd_C Pd_N Rb Rh Rh_N Ru_N Sc Sc_C Se_As Si_C Sn Ta_N Ti Ti_C Ti_N V V_N W W_C W_N Zn Zn_S Zn_Se Zn_Te Zr Zr_N
# BaO Be Ca CaO CaS CaSe Cd Co CrC CrN CsF FeN GaN InSb IrC IrN K KBr LaC LaN Li Mg MgO MgS MnC MnN MnO Mo Na NaF Nb NbN Ni NiC OsC OsN PdC PdN Rb Rh RhN RuN Sc ScC SeAs SiC Sn TaN Ti TiC TiN V VN W WC WN Zn ZnS ZnSe ZnTe Zr ZrN

# Pt LiH InP MoN Al NaCl ZrC Ge AgCl NiAl PtC CoN Sr AlAs FeC HfN CoAl RbI PtN BN GaAs Fe RhC FeAl HfC LiI NiN RuC MoC Ru CoC InAs BaSe GaP CdO BP TaC Si Pd GaSb LiF VC Ag Cu Au Ir NbC CsI Os Ta C MnS ScN LiCl BaO Be Ca CaO CaS CaSe Cd Co CrC CrN CsF FeN GaN InSb IrC IrN K KBr LaC LaN Li Mg MgO MgS MnC MnN MnO Mo Na NaF Nb NbN Ni NiC OsC OsN PdC PdN Rb Rh RhN RuN Sc ScC SeAs SiC Sn TaN Ti TiC TiN V VN W WC WN Zn ZnS ZnSe ZnTe Zr ZrN
# Pt Li_H In_P Mo_N Al Na_Cl Zr_C Ge Ag_Cl Ni_Al Pt_C Co_N Sr Al_As Fe_C Hf_N Co_Al Rb_I Pt_N B_N Ga_As Fe Rh_C Fe_Al Hf_C Li_I Ni_N Ru_C Mo_C Ru Co_C In_As Ba_Se Ga_P Cd_O B_P Ta_C Si Pd Ga_Sb Li_F V_C Ag Cu Au Ir Nb_C Cs_I Os Ta C Mn_S Sc_N Li_Cl
