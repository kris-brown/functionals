from typing import Tuple as T


def primcell(pth: str, struct: str, cellvol: float, volprimrat: float
             ) -> T[float, float, float]:
    ''' '''
    import os
    from ase.io import read
    import ase
    import scipy.optimize as opt

    # Compute volume to lattice constant^3 ratio
    sqrt2, sqrt3 = 2 ** (1 / 2), 3 ** (1 / 2)
    a_prim_corr = dict(fcc=sqrt2, bcc=2 * sqrt3 / 3, rocksalt=sqrt2,
                       diamond=sqrt2, zincblende=sqrt2)
    # a_conv_corr = dict(hcp=1, bcc=1, )
    with open(os.path.join(pth, 'latopt/OUTCAR'), 'r') as f:
        assert 'General timing' in f.read()
    atoms = read(os.path.join(pth, 'latopt/OUTCAR'))

    r = opt.root_scalar(f=lambda x: ase.atoms.Cell(
        atoms.get_cell() * x).volume - float(cellvol), x0=0.5, x1=2)
    atoms.set_cell(
        atoms.get_cell() * r.root, scale_atoms=True)

    a, _, _, alpha, _, _ = atoms.get_cell_lengths_and_angles()
    prim = bool(round(alpha) != 90)

    conventional_vol = float(cellvol) * (
        volprimrat if prim else 1)

    # if not prim and struct not in a_conv_corr:
    #     print(pth, struct, a)
    #     breakpoint()
    lattice_constant = a * (
        a_prim_corr[struct] if prim else 1)  # a_conv_corr[struct])

    volrat = conventional_vol / lattice_constant**3
    return volrat, conventional_vol, prim
