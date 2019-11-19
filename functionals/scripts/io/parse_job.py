from typing import List as L, Tuple as T, Set as S
from json import load,dumps

def parse_job(root:str, valid:S[str]
              )->T[L[str],L[int],L[int],L[float],L[str],L[str],
                   L[float],L[bool],L[float]]:
    from os      import listdir, environ, walk
    from os.path import join, exists
    from re      import findall,compile,MULTILINE, search
    from functionals.scripts.load.parse_incar import parse_incar
    from functionals.scripts.load.parse_contribs_vasp import parse_contribs_vasp
    from functionals.scripts.load.parse_eigenval import parse_eigenval

    bf = '/Users/ksb/functionals/data/beefs/beef.json'
    pths,pws,ecs,fxs,contribs,engs,intoccs,mags,rts = [],[],[],[],[],[],[],[],[] # type: ignore
    for pth, subdirList, fileList in walk(root):
        if 'INCAR' in fileList and 'bkup' not in pth: #any([search(v,pth) for v in valid]):
            incar = parse_incar(join(pth,'INCAR'))
            pths.append(pth)
            pws.append(int(incar['encut']));
            ecs.append(float(incar['ediff']));

            # Verify job completed successfully
            if exists(join(pth,'OUTCAR')) and ecs[-1] < 1e-4:
                with open(join(pth,'OUTCAR'),'r')  as f: outcar = f.read()
                with open(join(pth,'OSZICAR'),'r') as f: nsteps = sum(1 for _ in f)
            else: outcar =  ''

            if not ('General timing' in outcar and nsteps < 800): # completed
                engs.append(None); intoccs.append(None); contribs.append(None)
                rts.append(None); mags.append(None); fxs.append(None)
            else: # need to fill out eng/intocc/contrib/rt/mag
                if incar['metagga'] == 'BF':
                    # a1msb = [incar[x] for x in ['a1'+y for y in '12345']+['msb']]
                    with open(bf, 'r') as f: beef = load(f)
                    fxs.append(dumps(beef)) # we used to also store a1msb here
                    contribs.append(parse_contribs_vasp(outcar))
                else:
                    contribs.append(None) # type: ignore
                    if incar['metagga'] == 'SCAN':
                        fxs.append("SCAN")
                    elif incar['metagga'] is None and incar['gga']=='PE':
                        fxs.append('PBE')
                    else:
                        import pdb;pdb.set_trace(); raise ValueError()

                # GET ENERGY
                pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
                match = findall(pat, outcar); assert match
                engs.append(float(match[-1]))
                # Get runtime
                patt = r'Elapsed time \(sec\):\s+(\d+\.\d+)'
                match = findall(patt, outcar); assert match
                rts.append(int(float(match[-1])/60.))

                # GET OCCUPATION NUMBERS
                with open(join(pth,'EIGENVAL'),'r') as g:
                    io,mag = parse_eigenval(g.read())
                    intoccs.append(io); mags.append(mag)

    return pths,rts,pws,ecs,fxs,contribs,engs,intoccs,mags
