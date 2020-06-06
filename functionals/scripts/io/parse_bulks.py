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
    def err(x: str, engs: str = None, vols: str = None) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], str]:
        return None, None, None, vols, engs, None, x

    def read(car: str) -> Any:
        try:
            return io.read(car)
        except Exception as e:
            return "Exception:" + str(e)

    def process(pth: str) -> T[
            O[float], O[str],
            O[float], O[str], O[str], O[str], str]:

        latopt = read(os.path.join(pth, 'latopt/OUTCAR'))
        if latopt[:10] == "Exception:":
            return err("latopt - " + latopt)

        rpth = os.path.join(pth, 'latopt/run.log')
        if not os.path.exists(rpth):
            return err("No latopt/run.log")
        else:
            with open(rpth, 'r') as f:
                if 'fatal error' in f.read():
                    return err('fatal error in latopt run.log')

        # Analyze eos
        eos1engs, eos1vols = [], []
        for i in [-2, -1, 0, 1, 2]:
            strain = os.path.join(pth, 'eos/strain_%d' % i)
            if not os.path.exists(strain):
                return err("Missing  " + strain)
            _, [c_c_i], [c_e_i], [c_m_i], [z] = parse_job(strain)
            if z:
                return err(strain + ' ' + z)
            if i == 0:
                c_c = c_c_i
                c_e = c_e_i
                c_m = c_m_i
            if c_e_i is None:
                return err("Bad EOS1 job " + strain)
            else:
                eos1engs.append(c_e_i)
                eos1vols.append(io.read(strain + '/POSCAR').get_volume())
        for i, j in [(0, 1), (1, 2), (3, 2), (4, 3)]:
            if eos1engs[i] < eos1engs[j]:
                estr = "Bad EOS1 eng ordering %i %s - %i %s"
                return err(estr % (i, eos1engs[i], j, eos1engs[j]),
                           str(eos1engs), str(eos1vols))

        # Analyze eos2
        if not os.path.exists(os.path.join(pth, 'eos2')):
            vols = engs = conts = None
        else:
            vols_, engs_, conts_ = [], [], []
            for i in [-2, -1, 0, 1, 2]:
                strain = os.path.join(pth, 'eos2/strain_%d' % i)
                if not os.path.exists(strain):
                    return err("Missing  "+strain)
                _, [cont], [eng_], [y], [z] = parse_job(strain)
                if z:
                    return err(strain+' '+z)
                if eng_ is None:
                    return err("Bad EOS1 job "+strain)
                else:
                    engs_.append(eng_)
                conts_.append(json.loads(cont) if cont else [])
                vols_.append(io.read(strain+'/POSCAR').get_volume())

            if len(set(np.round(np.diff(vols_), 6))) > 1:
                return err('Irregular spacing of EOS %s' % np.diff(vols_))

            assert len(vols_) == 5
            vols, engs, conts = map(   # type: ignore
                json.dumps, [vols_, engs_, conts_])

            for i, j in [(0, 1), (1, 2), (3, 2), (4, 3)]:
                if engs_[i] < engs_[j]:
                    estr = "Bad EOS2 eng ordering %i %s - %i %s"
                    return err(estr % (i, engs_[i], j, engs_[j]))

        if isinstance(conts, str) and conts[:3] == '[[]':
            conts = None

        # Analyze phonon
        # TODO

        return c_e, c_c, c_m, vols, engs, conts, ''

    # Initialize
    pths, mats, xcs, engs, conts, mags, eosvols, eosengs, \
        eosconts, errs = [], [], [], [], [], [], [], [], [], []

    # Main loop
    for xc in os.listdir(root):
        for mat in os.listdir(os.path.join(root, xc)):
            if 'CdS' in mat or 'CdT' in mat:
                continue
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
