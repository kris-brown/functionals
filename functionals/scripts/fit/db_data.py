from typing import Set as S, List as L, Tuple as T
from json import dumps,loads

def db_data(a_ce:str,a_bm:str,a_lc:str,b_ce:str,b_bm:str,b_lc:str,ce_:str,bm_:str,lc_:str,name:str,ce_calc_:str) -> str:
    '''
    Assemble datasets given a database connection
    '''
    from functionals.fit.data import Datum,Data

    # Different datasets for each decay parameter value of VASP BEEF output
    datas = [set(),set(),set(),set(),set()] # type: L[S[Datum]]

    # '$' delimited arrays on a per-material basis
    arrs  = [a_ce,a_bm,a_lc,b_ce,b_bm,b_lc,name,ce_,bm_,lc_,ce_calc_]

    for ace_,abm_,alc_,bce_,bbm_,blc_,n,_ce,_bm,_lc,_ce_calc in zip(*map(lambda x: x.split('$'),arrs)):
        # Convert strings to arrays (of length 5) or floats if they are defined
        ace,abm,alc,bce,bbm,blc = map(lambda x: loads(x) if x else '',  [ace_,abm_,alc_,bce_,bbm_,blc_])
        ce, bm, lc, ce_calc     = map(lambda x: float(x) if x else None,[_ce,_bm,_lc,_ce_calc])
        # For each decay value:
        for i in range(5):
            if ace and ce and ce_calc:
                datas[i].add(Datum(n, 'ce', ace[i], bce[i], float(ce)))
            if abm and bm:
                datas[i].add(Datum(n, 'bm', abm[i], bbm[i], float(bm)))
            if alc and lc:
                datas[i].add(Datum(n, 'lc', alc[i], blc[i], float(lc)))

    return dumps([Data(d).to_list() for d in datas])
