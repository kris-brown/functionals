from typing import Tuple as T

def r2_fit(pth:str) -> T[float,float,float]:
    from functionals.fit.fit import FitResult
    fit = FitResult.from_pth(pth).full
    x = fit.xs[fit.opt]
    data = fit.train
    a,b,c = [data.r2(x,y) for y in ['ce','bm','lc']]
    return a,b,c
