from typing import FrozenSet as S,  Tuple as T, Optional as O
from json import loads, dumps
import numpy as np # type: ignore

def get_minjob(eos_bm : float, eng_ : str, vol_ : str, cont_ : str)->T[str,str,O[str],float,float]:
    '''
    Picks finite subset of all jobs to reproduce ASE eos bulk modulus
    '''
    bm = float(eos_bm)
    # conversion constant eV/A^6 to GPa
    bm_conv  = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    # deserialize and resort lists according to increasing volume
    eng,vol,cont = map(loads,[eng_,vol_,cont_])
    if not cont: cont = [[]]*len(eng_)
    else: assert len(eng)==len(vol) and len(eng)==len(cont)
    m = eng.index(min(eng))
    valid = set(range(len(vol)))
    close,mid,long = [frozenset(filter(lambda x: x in valid,[m+x*i for i in range(-2,3)])) for x in range(1,4)]
    sets = filter(lambda x: len(x)>4,[close, mid, long, close | mid,  mid|long,
                                      close | mid | long])

    def bmerr(inds:S[int])->T[float,float,float]:
        '''Fit bm to subset of data, report abs diff from EOSbm'''
        e,v      = [[x[i] for i in inds] for x in [eng,vol]]
        fit      = np.polyfit(v,e,2)
        bmguess  = -fit[1] * bm_conv # APPROXIMATE
        volguess = -fit[1]/2*(fit[0])
        return abs(bmguess - bm), bmguess, volguess

    errs = {s : bmerr(s) for s in sets}
    bestinds,(besterr,best_bmguess,best_volguess) = min(errs.items(),key=lambda x:x[1])
    inds = sorted(bestinds,key=lambda i:eng[i])
    oute,outv,outc = [dumps([x[i] for i in inds]) for x in [eng,vol,cont]]
    if not any(cont): outc = None # type: ignore
    return oute,outv,outc,best_bmguess,best_volguess
