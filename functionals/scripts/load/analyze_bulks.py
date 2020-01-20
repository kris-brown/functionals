from typing import Tuple as T, Optional as O
import numpy as np


def analyze_bulks(pth: str, vols: str, engs: str, lats: str, contribs: str,
                  mags: str, strains: str
                  ) -> T[str, int, int, str, str, int,
                         bool, float, O[float], O[float], float]:
    from json import loads
    from os.path import join
    from ase.io import read
    from ase.eos import EquationOfState
    from ase.units import kJ
    import warnings
    warnings.filterwarnings("ignore")

    s_vols, s_engs, s_lats, s_contribs, s_mags, s_strains = map(
        loads, [vols, engs, lats, contribs, mags, strains])
    has_contribs = len(s_contribs) == len(s_engs)
    zip_contribs = s_contribs if has_contribs else ['' for _ in s_engs]

    # Analysis of vectors
    # --------------------
    eos = EquationOfState(s_vols, s_engs)

    try:
        eosvol, _, eosbm = eos.fit()  # type: T[float,float,float]
    except Exception:
        eosvol = eosbm = 0.

    # Determine if minimum is centered
    tups = list(sorted(zip(s_vols, s_engs, zip_contribs, s_mags, s_strains))
                ) or [(0., 0., '', 0., 0)]
    minind = tups.index(min(tups, key=lambda x: x[1]))

    left = (minind < 2)
    right = (len(s_engs) - minind) < 3

    morejobs = -1 if left else (1 if right else 0)
    centered = not (left or right)

    gap = np.max(np.diff(sorted(s_strains))) if len(s_vols) > 1 else 0
    success = all([len(s_engs) > 5, centered,
                   gap < 2, eosvol > 0])

    # Things to do for an arbitrary strain
    # ------------------------------------
    atoms = read(join(pth, 'POSCAR'))
    n_atoms = len(atoms)
    elems = list(atoms.get_atomic_numbers())
    n_elems = len(set(elems))
    composition = str({elem: elems.count(elem) for elem in sorted(set(elems))})
    elemstr = ',%s,' % ','.join(map(str, sorted(elems)))
    eosbms = eosbm / kJ * 1.0e24 if success else None
    sd = pth[:pth.rfind('/')]
    # lattice parameter / V^-1/3
    latratio = np.linalg.norm(atoms.get_cell()[0])/atoms.get_volume()**(1/3)
    lat = latratio * eosvol**(1/3) if success else None
    # if n_elems==1 and 32 in elems: import pdb;pdb.set_trace()
    return (sd, n_atoms, n_elems,
            composition, elemstr, morejobs,
            success, eosbms, lat, eosvol if success else None, latratio)
