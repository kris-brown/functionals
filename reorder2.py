from typing import List
import sys
import os
import shutil
import ase.io
import scipy.optimize as opt
import numpy as np
import functionals.templates as temp

# Constants
j = os.path.join
s = 'strain_%d'
s2 = 'dirPOSCAR_0%d'
xs = ['CONTCAR', 'OUTCAR', 'OSZICAR', 'PROCAR', 'IBZKPT', 'CHG', 'CHGCAR',
      'DOSCAR', 'REPORT', 'vasprun.xml', 'run.log', 'EIGENVAL', 'XDATCAR',
      'PCDAT']
gt = temp.jinja_env.get_template


def main(pth: str, vol: float, space: float) -> None:

    def set_vol(strain: int, setvol: float) -> None:
        poscar = j(pth, s % strain, "POSCAR")
        atoms = ase.io.read(poscar)
        r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
            atoms.get_cell() * x).volume - setvol, x0=0.5, x1=2)
        atoms.set_cell(
            atoms.get_cell() * r.root, scale_atoms=True)
        ase.io.write(poscar, atoms)

    def vols(msg: str) -> List[float]:
        vs = [(round(ase.io.read(j(pth, s % i, "POSCAR")).get_volume(), 2), i)
              for i in range(-2, 3)]
        print(msg, vs)
        diffs = set(np.around(np.diff([x for (x, _) in vs]).tolist(), 3))
        assert len(diffs) == 1, "unequal spacing"
        return [float(a) for a, _ in vs]

    centr = j(pth, s % 0)
    if os.path.exists(j(pth, 'bkup')):
        shutil.rmtree(j(pth, 'bkup'))
    shutil.copytree(centr, j(pth, 'bkup'))
    # orig_vols = vols('before')

    set_vol(0, vol)

    for x in xs:
        try:
            os.remove(j(centr, x))
        except Exception:
            pass

    dirs = [centr]
    for i in [-2, -1, 1, 2]:
        ipth = j(pth, s % i)
        shutil.rmtree(ipth)
        shutil.copytree(centr, ipth)
        set_vol(i, vol + space * i)
        dirs.append(ipth)

    vols('after')

    bash = gt('subAll.jinja').render(gam='', dirs=dirs)
    bashpth = os.path.join(pth, 'subVASP.sh')
    with open(bashpth, 'w') as g:
        g.write(bash)

    os.system("""
    cd %s;
    pwd;
    chmod 755 %s;
    bsub -q suncat -W 20:00 -n 8 %s;
    cd -;
    """ % (pth, bashpth, bashpth))


def main2(pth: str, vol: float, space: float) -> None:

    def set_vol(strain: int, setvol: float) -> None:
        poscar = j(pth, s2 % strain, "POSCAR")
        atoms = ase.io.read(poscar)
        r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
            atoms.get_cell() * x).volume - setvol, x0=0.5, x1=2)
        atoms.set_cell(
            atoms.get_cell() * r.root, scale_atoms=True)
        ase.io.write(poscar, atoms)

    def vols(msg: str) -> List[float]:
        vs = [(round(ase.io.read(j(pth, s2 % i, "POSCAR")).get_volume(), 2), i)
              for i in range(5)]
        print(msg, vs)
        diffs = np.around(np.diff([x for (x, _) in vs]).tolist(), 3)
        assert len(set(diffs)) == 1, "unequal spacing" + str(diffs)
        return [float(a) for a, _ in vs]

    centr = j(pth, s2 % 2)
    # vols('before')

    set_vol(2, vol)

    for x in xs:
        try:
            os.remove(j(centr, x))
        except Exception:
            pass

    dirs = [centr]

    bfcar = '/nfs/slac/g/suncatfs/ksb/re42/re42_gas/propane/BEEFCAR'
    for i in [-2, -1, 1, 2]:
        ipth = j(pth, s2 % (i + 2))
        shutil.rmtree(ipth)
        shutil.copytree(centr, ipth)
        set_vol(i + 2, vol + space * i)
        dirs.append(ipth)

    for i in range(5):
        shutil.copyfile(bfcar, j(j(pth, s2 % i), 'BEEFCAR'))

    vols('after')

    bash = gt('subAll.jinja').render(gam='', dirs=dirs)
    bashpth = os.path.join(pth, 'subVASP.sh')
    with open(bashpth, 'w') as g:
        g.write(bash)

    os.system("""
    cd %s;
    pwd;
    chmod 755 %s;
    bsub -q suncat3 -W 30:00 -n 16 %s;
    cd -;
    """ % (pth, bashpth, bashpth))


if __name__ == '__main__':
    center = float(sys.argv[2])
    space = float(sys.argv[3])

    if False:
        assert os.path.exists(sys.argv[1] + '/strain_0/POSCAR')
        main(sys.argv[1], center, space)
    else:
        assert os.path.exists(sys.argv[1])
        main2(sys.argv[1], center, space)
