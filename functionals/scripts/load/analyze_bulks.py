from typing import Tuple as T, Optional as O


def analyze_bulks(pth: str
                  ) -> T[int, int, str, str, O[float], O[float], O[float]]:
    import os
    from ase.io import read

    atoms = read(os.path.join(pth, 'latopt/POSCAR'))
    n_atoms = len(atoms)
    elems = list(atoms.get_atomic_numbers())
    n_elems = len(set(elems))
    composition = str({elem: elems.count(elem) for elem in sorted(set(elems))})
    elemstr = ',%s,' % ','.join(map(str, sorted(elems)))
    try:
        with open(os.path.join(pth, 'latopt/OUTCAR'), 'r') as f:
            assert 'General timing' in f.read()
        atoms = read(os.path.join(pth, 'latopt/OUTCAR'))
        a = atoms.get_cell_lengths_and_angles()[0]
        vol = atoms.get_volume()
        volrat = vol/a**3
    except Exception:
        a = vol = volrat = None
    return (n_atoms, n_elems,  composition, elemstr, a,
            vol, volrat)
