from typing import List as L, Tuple as T,  Optional as O


def parse_job(root: str,
              ) -> T[L[str], L[O[str]], L[O[float]], L[O[float]], L[str]]:
    from os import walk
    from os.path import join, exists
    from re import findall, search
    from functionals.scripts.load.parse_incar import parse_incar
    from functionals.scripts.load.parse_contribs_vasp \
        import parse_contribs_vasp
    from functionals.scripts.load.parse_eigenval import parse_eigenval

    def err(x: str) -> T[O[str], O[float], O[float], str]:
        return None, None, None, x

    def process(pth: str) -> T[O[str], O[float], O[float], str]:
        incar = parse_incar(join(pth, 'INCAR'))
        if float(incar['ediff']) > 1e-5:
            return err('Ediff %s' % incar['ediff'])
        # Verify job completed successfully
        if all([exists(join(pth, 'OUTCAR')), exists(join(pth, 'OSZICAR'))]):
            with open(join(pth, 'OUTCAR'), 'r') as f:
                outcar = f.read()
            with open(join(pth, 'OSZICAR'), 'r') as f:
                nsteps = sum(1 for _ in f)
        else:
            return err('Job not started')

        if 'General timing' not in outcar:
            return err('Unconverged outcar')
        elif nsteps > 800:
            return err('Too many steps in OSZICAR %d' % nsteps)

        # Get contribs
        if incar['metagga'] == 'BF':
            try:
                cont = parse_contribs_vasp(outcar) or None
            except (IndexError, ValueError) as e:
                return err(str(e))
        else:
            cont = None

        # GET mag
        with open(join(pth, 'EIGENVAL'), 'r') as g:
            io, mag = parse_eigenval(g.read())
        if '/atoms/' in pth and not io:
            return err('Noninteger occupation numbers')
        if mag and '/atoms/' not in pth:
            pat = r'magnetization\s+([-]?\d+\.\d+)'
            pat2 = r'NIONS\s+=\s+(\d+)'
            natom = int(search(pat2, outcar).groups()  # type: ignore
                        [0])
            mags_ = list(findall(pat, outcar))
            magg = float(mags_[-2])/natom or None
        else:
            magg = None

        # GET ENERGY
        pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
        match = findall(pat, outcar)
        if not match:
            return err("Cannot parse energy %s - %s" % (match, pth))
        eng = float(match[-1])

        return cont, eng, magg, ''

    pths, contribs, engs, mags, errs = [
    ], [], [], [], []
    for pth, subdirList, fileList in walk(root):
        if 'INCAR' in fileList and 'bkup' not in pth:
            pths.append(pth)
            a, b, c, d = process(pth)
            contribs.append(a)
            engs.append(b)
            mags.append(c)
            errs.append(d)

    return pths, contribs, engs, mags, errs
