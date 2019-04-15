from typing import Tuple as T, List as L, Dict as D

def parse_fitparams(pth:str) -> T[L[str],L[float],L[float],L[float]]:
    from os.path import join
    from os      import environ
    from csv     import DictReader
    fitcols = ['consts','reg','bm_weight','lat_weight']
    vals = {c:[] for c in fitcols} # type: D[str,L[str]]
    types = dict(consts=str,reg=float,bm_weight=float,lat_weight=float)
    with open(pth,'r') as fi:
        reader = DictReader(fi)
        for row in reader:
            for k,v in vals.items(): v.append(types[k](row[k]))
    a,b,c,d = tuple([list(x) for x in vals.values()])
    return a,b,c,d # type: ignore
