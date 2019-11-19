from typing import Tuple as T

def mk_plt(name:str, pth : str, ce : float, bm : float, lc : float
          ) -> T[bool, float,float,float]:
    from functionals.fit.fit import FitResult
    from os import makedirs
    from os.path import join
    f = FitResult.from_pth(pth)
    m = dict(ce=float(ce),bm=float(bm),lc=float(lc))
    met   = {k:(1/v) if v!=0. else 0. for k,v in m.items()}
    cv    = f.plot_cv(met=met,show=False)
    trans = f.plot_transfer(met=met,show=False) if f.cv else ''
    res   = f.full.resid(step=f.full.opt,show=False)
    fx    = f.full.fxviz(opt=f.full.opt,show=False)
    err_ce, err_bm, err_lat = f.full.opt_err()
    dir = '/Users/ksb/scp_tmp/fitplt/'+name
    for plt,pname in zip([cv,trans,res,fx],['cv','trans','res','fx']):
        if plt:
            makedirs(dir,exist_ok=True)
            with open(join(dir,pname+'.html'),'w') as fi: fi.write(plt)
    return True, err_ce, err_bm, err_lat
