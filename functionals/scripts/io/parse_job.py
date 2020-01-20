from typing import List as L, Tuple as T, Set as S
from json import load, dumps


def parse_job(root: str, valid: S[str]
              ) -> T[L[str], L[int], L[int], L[float], L[str], L[str],
                     L[float], L[bool], L[float]]:
    from os import walk
    from os.path import join, exists
    from re import findall, search
    from functionals.scripts.load.parse_incar import parse_incar
    from functionals.scripts.load.parse_contribs_vasp import parse_contribs_vasp
    from functionals.scripts.load.parse_eigenval import parse_eigenval

    bf = '/Users/ksb/functionals/data/beefs/beef.json'
    pths, pws, ecs, fxs, contribs, engs, intoccs, mags, rts = [
    ], [], [], [], [], [], [], [], []
    for pth, subdirList, fileList in walk(root):
        if 'INCAR' in fileList and 'bkup' not in pth:
            incar = parse_incar(join(pth, 'INCAR'))
            pths.append(pth)
            pws.append(int(incar['encut']))
            ecs.append(float(incar['ediff']))

            # Verify job completed successfully
            if all([exists(join(pth, 'OUTCAR')), exists(join(pth, 'OSZICAR')),
                    ecs[-1] < 1e-4]):
                with open(join(pth, 'OUTCAR'), 'r') as f:
                    outcar = f.read()
                with open(join(pth, 'OSZICAR'), 'r') as f:
                    nsteps = sum(1 for _ in f)
            else:
                outcar = ''

            if not ('General timing' in outcar and nsteps < 800):  # completed
                engs.append(None)
                intoccs.append(None)
                contribs.append(None)
                rts.append(None)
                mags.append(None)
                fxs.append(None)
            else:  # need to fill out eng/intocc/contrib/rt/mag
                if incar['metagga'] == 'BF':
                    if exists(join(pth, 'BEEFCAR')):
                        with open(join(pth, 'BEEFCAR'), 'r') as f:
                            fx = list(map(float, f.read().split()))
                            assert len(fx) == 66, fx
                            fxs.append(dumps(fx[:-2]))
                    else:
                        assert '/beef/' in pth
                        with open(bf, 'r') as f:
                            beef = load(f)
                        # we used to also store a1msb here
                        fxs.append(dumps(beef))
                    try:
                        contribs.append(parse_contribs_vasp(outcar))
                    except (IndexError, ValueError):
                        print('\n\n'+pth+'\n\n')
                        assert False
                else:
                    contribs.append(None)  # type: ignore
                    if incar['metagga'] == 'SCAN':
                        fxs.append("SCAN")
                    elif incar['metagga'] is None:
                        if incar['gga'] == 'PE':
                            fxs.append('PBE')
                        elif incar['gga'] == 'PS':
                            fxs.append('PBESOL')
                        else:
                            raise ValueError()
                    else:
                        raise ValueError()

                # GET ENERGY
                pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
                match = findall(pat, outcar)
                assert match
                engs.append(float(match[-1]))
                # Get runtime
                patt = r'Elapsed time \(sec\):\s+(\d+\.\d+)'
                match = findall(patt, outcar)
                assert match
                rts.append(int(float(match[-1])/60.))

                # GET OCCUPATION NUMBERS
                with open(join(pth, 'EIGENVAL'), 'r') as g:
                    io, mag = parse_eigenval(g.read())
                    intoccs.append(io)
                    if mag and '/atoms/' not in pth:
                        pat = r'magnetization\s+([-]?\d+\.\d+)'
                        pat2 = r'NIONS\s+=\s+(\d+)'
                        natom = int(search(pat2, outcar).groups()
                                    [0])  # type: ignore
                        mags_ = list(findall(pat, outcar))
                        mag = float(mags_[-2])/natom

                    mags.append(mag)
    return pths, rts, pws, ecs, fxs, contribs, engs, intoccs, mags
