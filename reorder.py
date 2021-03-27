from typing import List
import sys
import os
import ase.io
import scipy.optimize as opt
import numpy as np

# Constants
j = os.path.join
s = 'strain_%d'
strains = [-2, -1, 0, 1, 2]

xs = ['CONTCAR', 'OUTCAR', 'OSZICAR', 'PROCAR', 'IBZKPT', 'CHG', 'CHGCAR',
      'DOSCAR', 'REPORT', 'vasprun.xml', 'run.log', 'EIGENVAL', 'XDATCAR',
      'PCDAT']


def main(root: str, vol: float, space: float) -> None:
    pth = root + '/strain_%d/'
    pos = pth + 'POSCAR'
    newvols = [vol + space * i for i in strains]

    def set_vol(strain: int, setvol: float) -> None:
        poscar = pos % strain
        atoms = ase.io.read(poscar)
        r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
            atoms.get_cell() * x).volume - setvol, x0=0.5, x1=2)
        atoms.set_cell(
            atoms.get_cell() * r.root, scale_atoms=True)
        ase.io.write(poscar, atoms)

    def vols(msg: str) -> List[float]:
        vs = [(round(ase.io.read(pos % i).get_volume(), 2), i)
              for i in strains]
        print(msg, vs)
        diffs = np.around(np.diff([x for (x, _) in vs]).tolist(), 3)
        if len(set(diffs)) != 1:
            print("unequal spacing" + str(diffs))
            breakpoint()
        return [float(a) for a, _ in vs]

    for strain, newvol in zip(strains, newvols):
        for x in xs:
            if os.path.exists(pth % strain + x):
                os.remove(pth % strain + x)
        set_vol(strain, newvol)

    vols('after')

    os.system(""" cd %s; pwd; bsub -q suncat3 -W 30:00 -n 16 %s; cd -;
    """ % (root, root + 'subVASP.sh'))


if __name__ == '__main__':
    _, xc, mat, center, space = sys.argv
    root = '/nfs/slac/g/suncatfs/ksb/beefjobs/bulks/%s/%s/eos/' % (xc, mat)
    assert os.path.exists(root + 'strain_0/POSCAR')
    main(root, float(center), float(space))
