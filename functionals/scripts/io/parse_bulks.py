from typing import Any, List as L, Tuple as T, Optional as O
import numpy as np


def parse_bulks(root: str
                ) -> T[L[str], L[str], L[str], L[O[float]], L[O[str]],
                       L[O[float]], L[O[str]], L[O[str]], L[O[str]],
                       L[O[float]], L[str]]:
    import os
    import ase.io as io
    import ase.units
    import json
    from functionals.scripts.io.parse_job import parse_job
    from ase.eos import EquationOfState
    from functionals.fit.data import allmats  # materials we care about

    # Helper functions

    def err(x: str, engs: str = None, vols: str = None) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], O[float], str]:
        return None, None, None, vols, engs, None, None, x

    def read(car: str) -> Any:
        try:
            return io.read(car)
        except Exception as e:
            return "Exception:" + str(e)

    def process(pth: str) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], O[float], str]:

        # Analyze eos
        if not os.path.exists(os.path.join(pth, 'eos')):
            return err("No EOS ")
        else:
            vols_, engs_, conts_ = [], [], []
            for i in [-2, -1, 0, 1, 2]:
                strain = os.path.join(pth, 'eos/strain_%d' % i)
                if not os.path.exists(strain):
                    return err("Missing  " + strain)
                _, [cont], [eng_], [mag_], [z] = parse_job(strain)
                if i == 0:
                    c_c, c_e, c_m = cont, eng_, mag_
                if z:
                    return err(strain + ' ' + z)
                if eng_ is None:
                    return err("Bad EOS job " + strain)
                else:
                    engs_.append(eng_)
                conts_.append(json.loads(cont) if cont else [])
                vols_.append(io.read(strain + '/POSCAR').get_volume())

            if len(set(np.round(np.diff(vols_), 3))) > 1:
                return err('Irregular spacing of EOS %s' % np.diff(vols_))

            assert len(vols_) == 5
            vols, engs, conts = map(
                json.dumps, [vols_, engs_, conts_])

            for i, j in [(0, 1), (1, 2), (3, 2), (4, 3)]:
                if engs_[i] < engs_[j]:
                    estr = "Bad EOS eng ordering %i %s - %i %s"
                    return err(estr % (i, engs_[i], j, engs_[j]), engs, vols)

        if isinstance(conts, str) and conts[:3] == '[[]':
            conts = None  # type: ignore

        volumes, energies = zip(*sorted(zip(vols_, engs_)))
        eos = EquationOfState(volumes, energies)
        v0, E0, B = eos.fit()
        B *= 1e24 / ase.units.kJ  # GPa
        return c_e, c_c, c_m, vols, engs, conts, B, ''

    # Initialize
    pths, mats, xcs, engs, conts, mags, eosvols, eosbms, eosengs, \
        eosconts, errs = [], [], [], [], [], [], [], [], [], [], []

    # Main loop
    for xc in filter(lambda x: '.' not in x, os.listdir(root)):
        for mat in filter(lambda x: x.replace('_', '') in allmats,
                          os.listdir(os.path.join(root, xc))):
            pth = os.path.join(root, xc, mat)
            pths.append(pth)
            mats.append(mat.replace('_', ''))
            xcs.append(xc)
            a, b, c, d, e, f, g, h = process(pth)
            engs.append(a)
            conts.append(b)
            mags.append(c)
            eosvols.append(d)
            eosengs.append(e)
            eosconts.append(f)
            eosbms.append(g)
            errs.append(h)
    return (pths, mats, xcs, engs, conts, mags, eosvols,
            eosengs, eosconts, eosbms, errs)


if __name__ == '__main__':
    print(*parse_bulks('/Users/ksb/scp_tmp/vauto/bulks'), sep='\n\n\n')
