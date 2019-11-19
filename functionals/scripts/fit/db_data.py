from typing import Set as S, List as L, Tuple as T
from json import dumps,loads

def db_data(ab_ce:str,ab_bm:str,ab_lc:str,ce_:str,bm_:str,lc_:str,name:str,ce_calc_:str) -> str:
    '''
    Assemble datasets given a database connection
    '''
    from functionals.fit.data import Datum,Data

    # Different datasets for each decay parameter value of VASP BEEF output
    data = set() # type: S[Datum]

    # '$' delimited arrays on a per-material basis
    arrs  = [ab_ce,ab_bm,ab_lc,name,ce_,bm_,lc_,ce_calc_]

    for abce_,abbm_,ablc_,n,_ce,_bm,_lc,_ce_calc in zip(*map(lambda x: x.split('$'),arrs)):
        # Convert strings to arrays (of length 5) or floats if they are defined
        abce, abbm, ablc = map(lambda x: loads(x) if x else '',  [abce_,abbm_,ablc_])
        ce, bm, lc, ce_calc     = map(lambda x: float(x) if x else None,[_ce,_bm,_lc,_ce_calc])
        if abce and ce and ce_calc:
            data.add(Datum(n, 'ce', abce[0], abce[1], float(ce)))
        if abbm and bm:
            data.add(Datum(n, 'bm', abbm[0], abbm[1], float(bm)))
        if ablc and lc:
            data.add(Datum(n, 'lc', ablc[0], ablc[1], float(lc)))

    return dumps(Data(data).to_list())
