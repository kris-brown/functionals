from typing import List as L, Tuple as T
from json import load,dumps

def parse_bulks(root:str)->T[L[str],L[str],L[int],L[int],L[float],L[str],
                             L[str],L[str],L[str],L[int],L[int],
                             L[str],L[str],L[float],L[float],
                             L[int],L[int],L[str],L[bool]]:

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

    try:

        # INITIALIZE VARS
        funroot  = environ['FUNCTIONALS_ROOT']

        pths,names,rts,pws,econvs,fxs, \
        contribs,engs,vols,n_atoms,n_elems,\
        comps,figs,eosbms,lats,\
        sls,shs,incs,sucs \
            = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

        randroot = environ['HOME']
        suffix   = 'tmp_'+''.join(choices(ascii_lowercase,k=8))+'.png'
        randpth  = join(randroot,suffix)

        # Traverse directory hierarchy with 3 levels (overall,batch,element)
        stordirs = listdir(root)
        for sd in stordirs: # 'date_description'
            folder = join(root,sd)
            for atom in listdir(folder):
                pth = join(folder,atom)
                strains = listdir(pth)

                # NAME AND PATH
                pths.append(pth); names.append(atom)

                s_engs     = [] # complete vector of energies
                s_vols     = [] # complete vector of volumes
                s_contribs = [] # complete vector of xc contributions
                s_lats     = [] # complete vector of lattice constants
                s_coms     = set() # list of strains that have completed
                rt,lo,hi   = 0,10000,-10000
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

                all = set(range(lo,hi+1)) # all strains
                incs.append(' '.join(map(str,sorted(all - s_coms))) or None)
                sls.append(lo); shs.append(hi); rts.append(rt)

                if len(s_vols) < 6:
                    inds,eos_vol = [0],0
                else:

                    # Analysis of vectors
                    #--------------------
                    eos = EquationOfState(s_vols,s_engs) # type: ignore

                    try:
                        eosvol,eoseng,eosbm = eos.fit() # type: T[float,float,float]
                    except:
                        eosvol,eoseng,eosbm = 0.,0.,0.

                    minind = s_engs.index(min(s_engs))
                    centered = (minind > 2) and (len(s_engs) - minind > 3)

                    if len(inds) > 5 and centered and eosvol > 0: # our bulk is truly complete, start adding to results

                        sucs.append(True)

                        # Make plot to file, then convert to string and clean up
                        eos.plot(randpth)
                        encoded = b64encode(open(randpth, "rb").read()).decode("utf-8")
                        figs.append(encoded); remove(randpth); plt.clf()

                        # Isolate the best 5 jobs,Convert vectors into strings
                        min_inds = [s_engs.index(x) for x in sorted(s_engs)[:5]]
                        lats.append(s_lats[min_inds[0]])
                        m_vols,m_engs,m_contribs = [[x[i] for i in min_inds]
                                                    for x in [s_vols,s_engs,s_contribs]] # filtered
                        vols.append(dumps(m_vols));engs.append(dumps(m_engs));
                        contribs.append(dumps(m_contribs))

                        # Things to do for an arbitrary strain
                        #------------------------------------

                        n_atoms.append(len(atoms))
                        elems = list(atoms.get_atomic_numbers())
                        n_elems.append(len(set(elems)))
                        comps.append(str({elem:elems.count(elem) for elem in sorted(set(elems))}))
                        # GET INCAR PARAMETERS
                        incar = parse_incar(join(dir,'INCAR'))
                        pws.append(incar['encut']); econvs.append(incar['ediff'])
                        eosbms.append(eosbm / kJ * 1.0e24)
                        # GET XC-SPECIFIC INFO
                        if incar['metagga'] == 'BF':
                            a1msb = [incar[x] for x in ['a1'+y for y in '12345']+['msb']]
                            with open(join(funroot,'data/beef.json'),'r') as f: beef = load(f)
                            fxs.append(dumps([beef]+a1msb))
                        else:
                            if incar['metagga'] == 'SCAN+rVV10':
                                fxs.append("SCAN")
                            elif incar['metagga'] is None and incar['gga']=='PE':
                                fxs.append('PBE')
                            else:
                                import pdb;pdb.set_trace()
                                raise ValueError()
                if not (len(inds) > 5 and eosvol > 0):
                    sucs.append(False)
                    pws.append(None); econvs.append(None);fxs.append(None) # type: ignore
                    contribs.append(None);engs.append(None);vols.append(None) # type: ignore
                    n_atoms.append(None); n_elems.append(None) # type: ignore
                    comps.append(None); figs.append(None); eosbms.append(None) # type: ignore
                    lats.append(None)


        return (pths,names,rts,pws,econvs,fxs,contribs,engs,vols,n_atoms,n_elems, # type: ignore
                comps,figs,eosbms,lats,sls,shs,incs,sucs)

    except Exception as e:
        import traceback,pdb
        traceback.print_exc(); print(pth,e);pdb.set_trace(); assert False
