from typing import Any, List as L, Tuple as T, Optional as O
import numpy as np


def parse_bulks(root: str
                ) -> T[L[str], L[str], L[str], L[O[float]], L[O[str]],
                       L[O[float]], L[O[str]], L[O[str]], L[O[str]], L[str]]:
    import os
    import ase.io as io
    import json
    from functionals.scripts.io.parse_job import parse_job

    # Helper functions
    def err(x: str) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], str]:
        return None, None, None, None, None, None, x

    def read(car: str) -> Any:
        try:
            return io.read(car)
        except Exception as e:
            return str(e)

    def process(pth: str) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], str]:
        # Analyze Latopt
        [c_p], [c_c], [c_e], [c_m], [c_er] = parse_job(
            os.path.join(pth, 'latopt'))
        if c_er:
            return err('latopt '+c_er)
        atoms = read(os.path.join(pth, 'latopt/OUTCAR'))
        if isinstance(atoms, str):
            return err('latopt ' + atoms)

        # Analyze eos
        vols, engs, conts = [], [], []
        for i in [-2, -1, 1, 2]:
            strain = os.path.join(pth, 'eos/strain_%d' % i)
            if not os.path.exists(strain):
                break
            _, [cont], [eng_], [y], [z] = parse_job(strain)
            if z:
                return err(strain+' '+z)
            engs.append(eng_)
            conts.append(json.loads(cont) if cont else [])
            vols.append(io.read(strain+'/POSCAR').get_volume())
            if i == -1:
                vols.append(atoms.get_volume())
                engs.append(c_e)
                conts.append(c_c)

        if len(set(np.round(np.diff(vols), 6))) > 1:
            return err('Irregular spacing of EOS')

        if len(vols) != 5:
            vols = engs = conts = None  # type: ignore
        else:
            vols, engs, conts = map(  # type: ignore
                json.dumps, [vols, engs, conts])

        if isinstance(conts, str) and conts[:3] == '[[]':
            conts = None

        # Analyze phonon
        # TBD

        return c_e, c_c, c_m, vols, engs, conts, ''  # type: ignore

    # Initialize
    pths, mats, xcs, engs, conts, mags, eosvols, eosengs, \
        eosconts, errs = [], [], [], [], [], [], [], [], [], []

    # Main loop
    for xc in os.listdir(root):
        for mat in os.listdir(os.path.join(root, xc)):
            pth = os.path.join(root, xc, mat)
            pths.append(pth)
            mats.append(mat)
            xcs.append(xc)
            a, b, c, d, e, f, g = process(pth)
            engs.append(a)
            conts.append(b)
            mags.append(c)
            eosvols.append(d)
            eosengs.append(e)
            eosconts.append(f)
            errs.append(g)
    return (pths, mats, xcs, engs, conts, mags, eosvols,
            eosengs, eosconts, errs)


if __name__ == '__main__':
    print(*parse_bulks('/Users/ksb/scp_tmp/vauto/bulks'), sep='\n\n\n')
