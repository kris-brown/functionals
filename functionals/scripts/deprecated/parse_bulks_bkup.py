from typing import List as L, Tuple as T
from json import load,dumps

def parse_bulks(root:str)->T[L[str],L[str],L[int],L[int],L[str],
                             L[str],L[str],L[str],L[str],L[str],L[int],L[int],
                             L[str],L[float],L[float],
                             L[int],L[int],L[int],L[str],L[bool]]:

    from os        import listdir, environ, remove
    from os.path   import join, exists
    from re        import findall
    from ase.io    import read  # type: ignore
    from string    import ascii_lowercase
    from random    import choices

    from base64    import b64encode
    from ase.eos   import EquationOfState   # type: ignore
    from ase.units import kJ                # type: ignore

    from   numpy        import polyfit,poly1d,mean,std,array # type: ignore
    from   numpy.linalg import norm # type: ignore
    import warnings; warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt # type: ignore

    from functionals.scripts.load.parse_incar import parse_incar
    from functionals.scripts.load.parse_contribs_vasp import parse_contribs_vasp
    bf = '/Users/ksb/functionals/data/beefs/beef.json'

    pths,names,rts,pws,fxs, \
    alleng,allvol,allcontrib,n_atoms,n_elems,\
    comps,eosbms,lats,\
    sls,shs,mjs,incs,sucs \
        = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    randroot = environ['HOME']
    suffix   = 'tmp_'+''.join(choices(ascii_lowercase,k=8))+'.png'
    randpth  = join(randroot,suffix)

    # Traverse directory hierarchy with 3 levels (overall,batch,element)
    stordirs = listdir(root)
    for sd in stordirs: # xc
        if sd[:3]== 'old': continue
        folder = join(root,sd)
        for bulk in listdir(folder):
            pth = join(folder,bulk)
            strains = listdir(pth)

            # NAME AND PATH
            s_engs     = [] # complete vector of energies
            s_vols     = [] # complete vector of volumes
            s_contribs = [] # complete vector of xc contributions
            s_lats     = [] # complete vector of lattice constants
            s_coms     = set() # list of strains that have completed
            rt,lo,hi   = 0,10000,-10000
            atoms      = None
            # Things to do for every strain
            for strain in strains:
                dir       = join(pth,strain)
                strainval = int(strain.split('_')[-1])
                lo = min(lo,strainval); hi = max(hi,strainval)

                # Verify job completed successfully
                if exists(join(dir,'OUTCAR')):

                    with open(join(dir,'OUTCAR'),'r')  as f: outcar = f.read()
                    with open(join(dir,'OSZICAR'),'r') as f: nsteps = sum(1 for _ in f) # counts # of lines in file
                else: outcar=''

                if 'General timing' in outcar and nsteps < 800: # completed
                    atoms = read(join(dir,'POSCAR'))

                    s_vols.append(atoms.get_volume())
                    s_lats.append(norm(atoms.get_cell()[0]))
                    # GET ENERGY
                    pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
                    match = findall(pat, outcar); assert match
                    s_engs.append(float(match[-1]))
                    s_coms.add(strainval)

                    # GET CONTRIBS
                    try:    s_contribs.append(parse_contribs_vasp(outcar))
                    except: s_contribs.append('')

                    # Get runtime THIS WILL BREAK IF WE HAVE NO COMPLETED JOBS?
                    patt = r'Elapsed time \(sec\):\s+(\d+\.\d+)'
                    match = findall(patt, outcar);
                    if match: rt = int(float(match[-1])/60.)
            if atoms is None: continue
            all = set(range(lo,hi+1)) # all strains
            incs.append(' '.join(map(str,sorted(all - s_coms))) or None)
            sls.append(lo); shs.append(hi); rts.append(rt)
            pths.append(pth); names.append(bulk)


            # Analysis of vectors
            #--------------------
            eos = EquationOfState(s_vols,s_engs) # type: ignore


            try:
                eosvol,eoseng,eosbm = eos.fit() # type: T[float,float,float]
            except:
                eosvol,eoseng,eosbm = 0.,0.,0.

            pairs = list(sorted(zip(s_vols,s_engs))) or [(0,0)]
            minind = pairs.index(min(pairs,key=lambda x:x[1]))
            left = (minind < 2)
            right = (len(s_engs) - minind) < 3
            mjs.append(-1 if left else (1 if right else 0))
            centered = not (left or right)
            #print('\n\n\n',pth)
            #import pdb;pdb.set_trace()

            suc =  bool((len(s_engs) > 5) and centered and (eosvol > 0))
            sucs.append(suc)

            # Isolate the best 5 jobs,Convert vectors into strings
            min_inds = [s_engs.index(x) for x in sorted(s_engs)[:5]]
            lats.append(s_lats[min_inds[0]] if min_inds else None)
            #m_vols,m_engs,m_contribs = [[x[i] for i in min_inds]
            #                            for x in [s_vols,s_engs,s_contribs]] # filtered
            #vols.append(dumps(m_vols));engs.append(dumps(m_engs));
            #contribs.append(dumps(m_contribs))
            allvol.append(dumps(s_vols))
            alleng.append(dumps(s_engs))
            allcontrib.append(dumps(s_contribs))

            # Things to do for an arbitrary strain
            #------------------------------------

            n_atoms.append(len(atoms))
            elems = list(atoms.get_atomic_numbers())
            n_elems.append(len(set(elems)))
            comps.append(str({elem:elems.count(elem) for elem in sorted(set(elems))}))
            # GET INCAR PARAMETERS
            incar = parse_incar(join(dir,'INCAR'))
            pws.append(int(incar['encut']))
            eosbms.append(eosbm / kJ * 1.0e24)

            # GET XC-SPECIFIC INFO
            if incar['metagga'] == 'BF':
                a1msb = [incar[x] for x in ['a1'+y for y in '12345']+['msb']]
                with open(bf,'r') as f: beef = load(f)
                fxs.append(dumps([beef]+a1msb))
            else:
                if incar['metagga'] == 'SCAN':
                    fxs.append("SCAN")
                elif incar['metagga'] is None and incar['gga']=='PE':
                    fxs.append('PBE')
                else:
                    import pdb;pdb.set_trace()
                    raise ValueError()
    print(list(map(len, [pths,names,rts,pws,fxs,alleng,allvol,allcontrib,n_atoms,n_elems, # type: ignore
            comps,eosbms,lats,sls,shs,mjs,incs,sucs])))
    import pdb;pdb.set_trace()

    return (pths,names,rts,pws,fxs,alleng,allvol,allcontrib,n_atoms,n_elems, # type: ignore
            comps,eosbms,lats,sls,shs,mjs,incs,sucs)

    # except Exception as e:
    #     import traceback,pdb
    #     traceback.print_exc(); print(pth,e);pdb.set_trace(); assert False
