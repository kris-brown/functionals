from typing import Set, Dict, Any, Tuple, Optional
import argparse
import collections
import csv
import json
import math
import os
import os.path
import subprocess
import ase
import ase.data
import ase.io
import numpy as np
import functionals.templates as temp

########
# DATA #
########


def str2bool(v: str) -> bool:
    return v.lower() in ("yes", "true", "t", "1")


ppth = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/potpaw_PBE'
root = '/nfs/slac/g/suncatfs/ksb/beefjobs/'
exe_root = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/'\
           'olsherlock/vasp.5.4.1beefcar/bin/vasp_'

data_root = '/'+os.path.join(*__file__.split('/')[:-3])+'/data'
struct_root = data_root + '/structures'
vol_csv = data_root+'/opt_vol.csv'
beef_root = data_root + '/beefs'

beefs = [x[:-5] for x in os.listdir(beef_root)]

atommag, n_electrons = {}, {}
Mat = collections.namedtuple(
    'Mat', ['name', 'struct', 'lat', 'bm', 'mag', 'ca_rat'])

with open(data_root+'/elements.csv', 'r') as f:
    r = csv.reader(f)
    next(r)
    for e, n, m in r:
        atommag[e], n_electrons[e] = int(m), int(n)


with open(data_root + '/expt.csv', 'r') as f:
    r = csv.reader(f)
    has_data = {x: (bool(a), bool(b), bool(c), bool(d))
                for x, a, b, c, d in r}

with open(data_root+'/initialguess.csv', 'r') as f:
    r = csv.reader(f)
    next(r)
    matdata = {k: Mat(k, s, float(l), float(b), float(m or '0'),
                      float(ca or '0'))
               for k, s, l, b, m, ca in r
               if k in has_data}

xcdict = dict(pbe=dict(gga='PE'), scan=dict(metagga='SCAN'),
              pbesol=dict(gga='PS'),
              )  # type: Dict[str,Dict[str,Any]]
xcdict.update(**{beef: dict(metagga='BF', lbeefens=True)
                 for beef in beefs if beef not in xcdict})


def only_ce(m: Mat) -> bool:
    ce, bm, lat, _ = has_data[m.name]
    return ce and not (bm or lat)

# def opt_vol() -> Dict[Tuple[str, str], Tuple[float, bool, bool, bool]]:
#     with open(vol_csv, 'r') as f:
#         dic = {(a, b): (float(c), str2bool(x), str2bool(y), str2bool(z))
#                for a, b, c, _, x, y, z in csv.reader(f)}
#     return dic


exedic = collections.defaultdict(lambda: 'vasp')  # type: Dict[str, str]
exedic['scan'] = 'scan'
##############################################################################
# HELPER FUNCTIONS #
n_dict = dict(fcc=4, bcc=2, hcp=2, cesiumchloride=2, diamond=8,
              zincblende=8, rocksalt=8)


def mk_traj(m: Mat) -> ase.Atoms:
    '''Create an Atoms object from name, structure, and lattice constant.'''
    c = [[1., 0, 0], [0, 1, 0], [0, 0, 1]]

    if m.struct == 'fcc':
        n = 4
        p = [[0, 0, 0], [0, 1/2, 1/2], [1/2, 0, 1/2], [1/2, 1/2, 0]]
    elif m.struct in ['bcc', 'cesiumchloride']:
        n = 2 if m.struct == 'bcc' else 1
        p = [[0, 0, 0], [1/2, 1/2, 1/2]]
    elif m.struct == 'hcp':
        n = 2
        p = [[0, 0, 0], [1/3, 2/3, 1/2]]
        c = [[1, 0, 0], [-0.5, 3**0.5/2, 0], [0, 0, m.ca_rat]]  # DIFFERENT
    elif m.struct in ['diamond', 'zincblende']:
        n = 8 if m.struct == 'diamond' else 4
        p = [[0, 0, 0], [1/4, 1/4, 1/4], [0, 1/2, 1/2], [1/4, 3/4, 3/4],
             [1/2, 0, 1/2], [3/4, 1/4, 3/4], [1/2, 1/2, 0], [3/4, 3/4, 1/4]]
    elif m.struct == 'rocksalt':
        n = 4
        p = [[0, 0, 0], [1/2, 0, 0], [0, 1/2, 1/2], [1/2, 1/2, 1/2],
             [1/2, 0, 1/2], [0, 0, 1/2], [1/2, 1/2, 0], [0, 1/2, 0]]
    else:
        raise NotImplementedError(m.struct)
    return ase.Atoms(m.name*n, cell=np.array(c)*m.lat, scaled_positions=p)


def onezero(one: int, zero: int) -> str:
    return ('%d*1.0' % one if one else '') + ' %d*0.0' % zero  # FERWE/FERDO


def get_strain(mat: Mat, strain_int: int) -> float:
    '''Gives the *lattice* strain, in order to get a *volume* strain.'''
    spacing = float(0.25/(mat.bm)**(1./3.))
    if mat in ['K', 'Ca', 'Rb', 'Sr']:
        spacing*3
    return (1 + strain_int * spacing)**(1/3)


def current_jobs() -> Set[str]:
    cmd = 'bjobs -o "job_name" -noheader'
    return set(subprocess.check_output(cmd, encoding='UTF-8',
                                       shell=True).split())


def readfile(pth: str) -> str:
    """Convert filepath to string of contents."""
    with open(pth, 'r') as f:
        return f.read()


def done(spth: str) -> bool:
    '''Identify if a job has already successfully completed.'''
    ozpth, outpth = [os.path.join(spth, f) for f in ['OSZICAR', 'OUTCAR']]

    don = os.path.exists(ozpth) and (len(readfile(ozpth).split('\n')) < 800)\
        and ('General timing' in readfile(outpth))
    return don


def mk_incar(d: Dict[str, Any], pth: str) -> None:
    """Convert dictionary to INCAR file."""
    def fmt(val: Any) -> str:
        if isinstance(val, bool):
            return '.TRUE.' if val else '.FALSE.'
        else:
            return str(val)

    with open(pth, 'w') as f:
        lines = ['{} = {}'.format(k.upper(), fmt(v)) for k, v in d.items()]
        f.write('\n'.join(lines)+'\n')


def potpos(pth: str, atoms: ase.Atoms) -> None:
    poscar, potcar = [os.path.join(pth, x) for x in ['POSCAR', 'POTCAR']]
    ase.io.write(poscar, atoms)
    elems = subprocess.check_output('head -n 1 '+poscar, encoding='UTF-8',
                                    shell=True).split()
    potcmd = ('cat ' + ' '.join(['%s/%s/POTCAR' % (ppth, x) for x in elems])
              + ' > ' + potcar)
    os.system(potcmd)


def attempt(retry: bool, pth: str, curr: Set[str]) -> bool:
    notdone = not done(pth)
    notcurr = os.path.join(pth, 'subVASP.sh') not in curr
    return retry or (notdone and notcurr)


def magdic(mat: Mat) -> Dict[str, str]:
    if mat.mag:
        mm = ' '.join([str(mat.mag)]*n_dict[mat.struct])
        return {'ispin': "2", 'magmom': mm}
    else:
        return dict(ispin='1')

##############################################################################
# MAIN FUNCTIONS #


def submit_atom(name: str, time: int, xc: str, retry: bool,
                incar: Dict[str, Any],
                curr: Set[str], orbitals: Optional[bool], sunc: str,
                beef: Optional[str]) -> None:
    assert sunc == '3'
    elem = ase.data.chemical_symbols.index(name)
    atoms = ase.Atoms(numbers=[elem], positions=[[1., 2., 3.]],
                      cell=[[19., 2., 1.],  [2., 20, 1.],  [5., 1, 21.]])
    magmom, elec = atommag[name], n_electrons[name]
    nb = max(8, int(math.ceil(0.6*elec+magmom if magmom else elec/2+0.5*1)))
    ferwedo = orbitals if (orbitals is not None) else (elem > 56)
    if ferwedo and magmom:
        nup, ndn = (elec + magmom)//2, (elec - magmom)//2
        magdic = dict(ismear=-2, ferwe=onezero(nup, nb - nup), ispin=2,
                      ferdo=onezero(ndn, nb - ndn))

    else:
        magdic = dict(sigma=0.0001, nupdown=magmom, ispin=2 if magmom else 1)

    magcar = dict(
        ldiag=False, isym=-1, amix=0.1 if xc == 'scan' else 0.2,
        bmix=0.0001, nbands=nb, **magdic)

    # Submit
    pth = os.path.join(root, 'atoms', xc, name)
    os.makedirs(pth, exist_ok=True)
    if attempt(retry, pth, curr):
        setup(pth=pth, atoms=atoms, incar={**incar, **magcar},
              kpts=(1, 1, 1), beef=beef)
        bash = temp.jinja_env.get_template(
            'subVASP.jinja').render(exe=exedic[xc], scan=xc == 'scan')
        _sub(pth, bash, time, sunc)


def phonon(mat: Mat, time: int, xc: str,
           retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
           beef: Optional[str]) -> None:
    raise NotImplementedError


def latopt(mat: Mat, time: int, xc: str,
           retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
           beef: Optional[str]) -> None:
    new = dict(ediff=1e-5, nsw=500, ediffg=-0.005, addgrid=True,
               ismear=0, ivdw=0, lreal=False, ibrion=2, isif=7)
    atoms = mk_traj(mat)
    pth = os.path.join(root, 'bulks', xc, mat.name, 'latopt')
    os.makedirs(pth, exist_ok=True)
    if attempt(retry, pth, curr):
        setup(pth=pth, atoms=atoms, incar={**incar, **new, **magdic(mat)},
              kpts=(10, 10, 10), beef=beef)
        bash = temp.jinja_env.get_template('subVASP.jinja').render(
            exe=exedic[xc], scan=xc == 'scan')
        _sub(pth, bash, time, sunc)


def eos(mat: Mat, time: int, xc: str,
        retry: bool, incar: Dict[str, Any], curr: Set[str], sunc: str,
        beef: Optional[str]) -> None:
    if only_ce(mat):
        print('Shouldnt do EOS on %s without bm/lat data' % mat.name)
        return None
    matpth = os.path.join(root, 'bulks', xc, mat.name, 'eos')
    os.makedirs(matpth, exist_ok=True)
    outpth = matpth[:-3]+'latopt'
    if not done(outpth):
        print('Cannot do EOS for %s - incomplete OUTCAR' % matpth)
        return None
    try:
        orig_atoms = ase.io.read(outpth+'/OUTCAR')
    except Exception:
        print('Cannot do EOS for %s - incomplete OUTCAR' % matpth)
        return None
    # Bulk-specific INCAR settings
    magcar = dict(sigma=0.01, ismear=0, **magdic(mat))

    dirs = []

    for strain_i in [-2, -1, 1, 2]:
        # Strained Traj
        lat_strain = get_strain(mat, strain_i)
        atoms = orig_atoms.copy()
        atoms.set_cell(atoms.get_cell()*lat_strain, scale_atoms=True)
        # Determine whether to submit again or not
        pth = os.path.join(matpth, 'strain_%d' % strain_i)
        os.makedirs(pth, exist_ok=True)
        if not done(pth):
            setup(pth=pth, atoms=atoms, incar={**incar, **magcar},
                  kpts=(10, 10, 10), beef=beef)
            dirs.append(pth)

    if attempt(retry, matpth, curr) and curr:
        gt = temp.jinja_env.get_template
        bash = gt('subAll.jinja').render(exe=exedic[xc], dirs=dirs,
                                         rm=xc != 'pbe')
        _sub(matpth, bash, time, sunc)


def _sub(pth: str, bash: str, time: int, sunc: str) -> None:
    bashpth = os.path.join(pth, 'subVASP.sh')
    with open(bashpth, 'w') as g:
        g.write(bash)
    os.chdir(pth)
    os.system('chmod 755 '+bashpth)
    args = [16 if sunc else 8, time, sunc, bashpth]
    cmd = 'bsub -n {} -W{}:09 -q suncat{} {}'.format(*args)
    os.system(cmd)


def setup(pth: str, atoms: ase.Atoms, incar: Dict[str, Any],
          kpts: Tuple[int, int, int], beef: Optional[str]) -> None:
    '''General setup of directory for atomic or bulk (at a strain).'''
    # Copy PBE wavecar if possible
    # pbejob = pth.replace('/beef/', '/pbe/').replace('/scan/', '/pbe/')
    # if all(['pbe' not in pth, 'opt' not in pth, done(pbejob),
    #         not os.path.exists(os.path.join(pth, 'WAVECAR'))]):
    #     os.system('cp {}/WAVECAR {}/WAVECAR'.format(pbejob, pth))

    # Write INCAR
    mk_incar(incar, os.path.join(pth, 'INCAR'))

    # Write BEEFCAR
    if beef:
        with open(os.path.join(pth, 'BEEFCAR'), 'w') as f:
            f.write(beef)

    # KPOINTS
    with open(os.path.join(pth, 'KPOINTS'), 'w') as g:
        g.write(temp.jinja_env.get_template('KPOINTS.jinja').render(kpts=kpts))

    # Write POSCAR and POTCAR
    potpos(pth, atoms)


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
                    type=lambda x: x.split())

typedict = dict(latopt=latopt, eos=eos, phonon=phonon)  # type: Any


def main() -> None:
    '''Submit either bulk or element singlepoint calculations.'''
    # Get/validate command line args
    args = parser.parse_args()
    for xc in args.xc:
        if xc not in ['pbe', 'scan', 'pbesol']:
            assert xc in beefs
            with open(os.path.join(beef_root, xc+'.json'), 'r') as f:
                beef = ' '.join(map(str, json.load(f)+[0., 1.]))
        else:
            beef = None  # type: ignore
        assert args.time > 0
        curr = current_jobs()

        incar = dict(encut=900., ediff=1e-5, algo='A', istart=2,
                     addgrid=True, ibrion=-1, npar=1 if args.mat else 4,
                     lasph=True, nelm=800,
                     lcharg=True, lwave=True, prec='ACCURATE', **xcdict[xc])

        common = dict(incar=incar, time=args.time, xc=xc, curr=curr,
                      retry=args.retry, sunc=args.sunc, beef=beef)
        if args.mat:
            assert (not args.elems) and args.type, args
            for m in args.mat:
                typedict[args.type](mat=m, **common)
        else:
            assert args.elems
            for e in args.elems:
                submit_atom(name=e, orbitals=args.orb, **common)


if __name__ == '__main__':
    main()
