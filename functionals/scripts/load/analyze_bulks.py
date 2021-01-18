from typing import Tuple as T


def analyze_bulks(pth: str, mat: str
                  ) -> T[int, int, str, str, str]:
    import os
    from ase.io import read
    from functionals.CLI.submit import matdata

    # Get an arbitrary POSCAR if possible
    try:
        struct = matdata[mat].struct
        try:
            atoms = read(os.path.join(pth, 'eos2/strain_0/POSCAR'))
        except Exception:
            atoms = read(os.path.join(pth, 'latopt/POSCAR'))
    except Exception:
        print('\n\n\n\nWeird pth ', pth, '\n\n\n')
        return (0, 0, '', '', '')

    # Compute facts about elemental composition
    n_atoms = len(atoms)
    elems = list(atoms.get_atomic_numbers())
    n_elems = len(set(elems))
    composition = str({elem: elems.count(elem) for elem in sorted(set(elems))})
    elemstr = ',%s,' % ','.join(map(str, sorted(elems)))

    return (n_atoms, n_elems, composition, elemstr, struct)
