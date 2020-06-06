from typing import List, Callable
import argparse
import os
import ase
from ase.build import fcc111, hcp0001, add_adsorbate
from functionals.CLI.submit import (mk_incar, root, potpos, matdata,
                                    _sub,
                                    xcdict, mk_traj)
import functionals.templates as temp

###########################################################

fccs = ['Rh', 'Pd', 'Pt', 'Cu']
hcps = ['Co', 'Ru']
sites = ['fcc', 'ontop']
ca_rats = dict(Co=1.6231901, Ru=1.582283159)

co = ase.Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.2)])


common = dict(encut=700., sigma=0.05, prec='ACCURATE', nelm=400,
              ediff=1e-5, nsw=500, ediffg=-0.01, addgrid=True,
              ismear=0, ivdw=0, lreal='A', ibrion=2)

###########################################################


def mksurf(mat: str, xc: str, site: str = None) -> ase.Atoms:
    '''Make a surface, optionally with CO on a site.'''
    bulkd = os.path.join(root, 'surf', xc, mat, 'bulk/OUTCAR')
    bulkdim = ase.io.read(bulkd).get_cell_lengths_and_angles()
    a, c = bulkdim[0], bulkdim[2]
    kwargs = dict(size=(2, 2, 4), vacuum=15.0, a=a)
    if mat in fccs:
        slab = fcc111(mat, **kwargs)
    else:
        slab = hcp0001(mat, c=c, **kwargs)
    slab.set_constraint(ase.constraints.FixAtoms(indices=range(8)))
    if site:
        add_adsorbate(slab=slab, adsorbate=co, height=1.5, position=site)
        slab.center(vacuum=15.0, axis=2)
    return slab


def subbulk(mat: str, xc: str) -> None:
    pth = os.path.join(root, 'surf', xc, mat, 'bulk')
    os.makedirs(pth, exist_ok=True)
    traj = mk_traj(matdata[mat])
    incar = dict(isif=3, **common)
    mk_incar(incar, os.path.join(pth, 'INCAR'))
    potpos(pth, traj)
    with open(os.path.join(pth, 'KPOINTS'), 'w') as g:
        g.write(temp.jinja_env.get_template(
            'KPOINTS.jinja').render(kpts=(10, 10, 10)))
    bash = temp.jinja_env.get_template('subVASP.jinja'
                                       ).render(rm=True)
    _sub(pth, bash, 10, '')


def submit(mat: str, xc: str, site: str = None,
           time: int = None) -> None:
    traj = mksurf(mat=mat, xc=xc, site=site)
    incar = dict(
        isym=0, npar=4, **common,
        **xcdict[xc])
    if mat == 'Co':  # ONLY MAG SURF
        incar.update(ispin=2, magmom='%d*1.8' % len(traj))
    else:
        incar.update(ispin=1)
    pth = os.path.join(root, 'surf', xc, mat, site or 'bare')
    os.makedirs(pth, exist_ok=True)
    # INCAR
    mk_incar(incar, os.path.join(pth, 'INCAR'))
    # POSCAR AND POTCAR
    potpos(pth, traj)
    # KPOINTS
    with open(os.path.join(pth, 'KPOINTS'), 'w') as g:
        g.write(temp.jinja_env.get_template(
            'KPOINTS.jinja').render(kpts=(4, 4, 1)))

    bash = temp.jinja_env.get_template('subVASP.jinja'
                                       ).render(rm=True)
    _sub(pth, bash, time or 20, '3')


###########################################################
def parse(opts: List[str]) -> Callable[[str], List[str]]:
    def f(x: str) -> List[str]:
        if x == 'all':
            return opts
        else:
            return x.split()
    return f


parser = argparse.ArgumentParser(description='Submit some jobs',
                                 allow_abbrev=True)
parser.add_argument('--xc', help='Functional', default='',
                    type=parse(['pbe', 'rpbe', 'msurf']))
parser.add_argument('--mat', help='Material', default='',
                    type=parse(fccs + hcps))
parser.add_argument('--bulk', help='Put anything to do bulk', default='')
parser.add_argument('--time', default=10, type=int, help='Walltime')

###########################################################
if __name__ == '__main__':
    args = parser.parse_args()
    for mat in args.mat:
        for xc in args.xc:
            if args.bulk:
                subbulk(mat, xc)
            else:
                for site in sites:
                    submit(mat=mat, xc=xc, site=site, time=args.time)
