from typing import List as L, Tuple as T, Optional as O, Any
import numpy as np

def analyze_bulks(pth: str, vols:str, engs:str, lats:str, contribs: str,
                  mags: str,strains: str
                )->T[str,int,int,
                     str,str,int,
                     bool,float,float,
                     float,float,float, float, str,
                     float,float,str,
                     str, O[str]]:
    from json      import loads, dumps
    from os.path   import join
    from ase.io    import read  # type: ignore
    from ase.eos   import EquationOfState   # type: ignore
    from ase.units import kJ                # type: ignore
    import warnings; warnings.filterwarnings("ignore")


    s_vols, s_engs, s_lats, s_contribs, s_mags, s_strains = map(loads, [vols, engs, lats, contribs, mags, strains])
    has_contribs = len(s_contribs) == len(s_engs)
    zip_contribs = s_contribs if has_contribs else ['' for _ in s_engs]
    # conversion constant eV/A^6 to GPa
    ev_a3_to_gpa  = (10**-9) * (1.602 * 10**-19) * (10**10)**3


    # Analysis of vectors
    #--------------------
    eos = EquationOfState(s_vols,s_engs) # type: ignore

    try:
        eosvol,eoseng,eosbm = eos.fit() # type: T[float,float,float]
    except:
        eosvol,eoseng,eosbm = 0.,0.,0.

    # Determine if minimum is centered
    tups = list(sorted(zip(s_vols,s_engs,zip_contribs,s_mags,s_strains))) or [(0.,0.,'',0.,0)]
    minind = tups.index(min(tups, key=lambda x:x[1]))
    minstrain = s_strains[minind]
    mag = tups[minind][3]

    left = (minind < 2)
    right = (len(s_engs) - minind) < 3

    morejobs = -1 if left else (1 if right else 0)
    centered = not (left or right)
    contig = centered and set(s_strains[minind-2:minind+3])==set(range(minstrain-2,minstrain+3))
    success =  bool((len(s_engs) > 5) and contig and (eosvol > 0))
    if success:
        minjobs = [tups[i] for i in range(minind-2,minind+3)]
        fit_vols, fit_engs, fit_contribs, _, _ = map(list, zip(*minjobs))
        fit = np.polyfit(fit_vols,fit_engs, 2)
        vol = eosvol # float(-fit[1]/(2*fit[0]))
        dx = np.mean(np.diff(fit_vols))
        assert np.max(np.diff(fit_vols)-dx) < 0.0001, vols  # check even spacing
        stencil = (-1*fit_engs[4] + 16*fit_engs[3] - 30*fit_engs[2] + 16*fit_engs[1] - fit_engs[0])/(12*dx**2)  #type: ignore
        bulkmod = stencil * vol * ev_a3_to_gpa # APPROXIMATE
        bulkmod_lstsq  = float(-fit[1]) * ev_a3_to_gpa # APPROXIMATE
        eng = float(np.dot(fit,[vol**2, vol, 1]))
        latratio = float(s_vols[0]**(1/3)/s_lats[0])
        lat = latratio * vol**(1/3)
        str_fv, str_fe, str_fc, abc = map(dumps, [fit_vols, fit_engs, fit_contribs,fit.tolist()])

    else:
        str_fv=str_fe=str_fc=bulkmod=bulkmod_lstsq=vol=eng=lat=latratio=abc=None # type: ignore

    # Things to do for an arbitrary strain
    #------------------------------------
    atoms = read(join(pth,'POSCAR'))
    n_atoms = len(atoms)
    elems = list(atoms.get_atomic_numbers())
    n_elems = len(set(elems))
    composition = str({elem:elems.count(elem) for elem in sorted(set(elems))})
    elemstr = ',%s,'%','.join(map(str,sorted(elems)))
    eosbms = eosbm / kJ * 1.0e24
    sd = pth[:pth.rfind('/')]

    return (sd, n_atoms, n_elems,
            composition, elemstr, morejobs,
            success,  eng, lat,
            vol, latratio, bulkmod, bulkmod_lstsq, abc,
            eosbms, mag,  str_fe,
            str_fv, str_fc if has_contribs else None)
