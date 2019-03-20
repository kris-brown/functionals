from typing import Tuple as T, List as L, Optional as O
from math   import log10
from ast    import literal_eval
from os     import environ
from json   import loads
import numpy as np # type: ignore

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, COALESCE, COUNT, EQ,
                    flatten, Env, Import, defaultEnv, Const, Constraint, JPath)

from functionals.fit import Fit as FitObj , FromMatrix #Functional,BEEF

from functionals.scripts.load.parse_fitstep   import  parse_fitstep
from functionals.scripts.fit.fitting          import  fitting
from functionals.scripts.fit.h_norm_const     import  one_time_function

################################################################################
################################################################################
################################################################################

functionals_db = '/Users/ksb/Documents/JSON/functionals.json'

################################################################################


################################################################################
def sqldict(d : dict) -> Expr:
    '''Turn a Python dict into a SQL (serialized) dict'''
    args = [[Lit(','), Lit(' "%s" : ['%k), GC(v), Lit(']')] for k,v in d.items()] # type: L[L[Expr]]
    return C(Lit('{'),*flatten(args)[1:],Lit('}'))

################################################################################
################################################################################
################################################################################

def fit(mod : Model) -> None:
    # Extract tables
    tabs = ['fit','calc','expt','species','fitparams']
    Fit, Calc, Expt, Species, Fitparams = map(mod.get, tabs) # type: ignore

    fit__calc, fit__fitparams = map(mod.get_rel,[Fit.r('calc'),Fit.r('fitparams')])

    fbq = Query(dict(f = Fit.id(),
                     p = Fit['pth'](),
                     r = Fit['result'](),
                     a = Fit['a1'](),
                     m = Fit['msb']()))

    ############################################################################

    resid_cols = ['r2_ce','r2_bm','r2_lat','c_viol']

    f_env = defaultEnv + Env(Import('functionals.fit.fit',Fit='FitObj'))

    def f_resid(db:str,x:str,c:int,d:int)->T[float,float,float,float]:
        return FitObj.costs(db,x,c,d)

    fc = JPath("calc", [fit__calc])
    fp = JPath('fitparams',[fit__fitparams])

    rq = Query(exprs = dict(f   = Fit.id(),
                            x = Fit['result'](),
                            c = Calc.id(fc),
                            d = Fit['decay']()),
                basis = ['fit'])
    db = '/Users/ksb/Documents/JSON/functionals.json'
    rpb = PyBlock(f_resid,
                   env  = f_env,
                   args = [Const(db),rq['x'],rq['c'],rq['d']],
                   outnames = resid_cols)

    resid =                                                                     \
        Gen(name    = 'resid',
            desc    = 'Computes metrics for a fit result',
            query   = rq,
            tags    = ['fit','parallel','DEP expt.coefs','DEP expt.expt_volume'],
            funcs   = [rpb],
            actions = [Fit(fit=rq['f'],
                           **{x:rpb[x] for x in resid_cols})])

    ########################################################################

    def f_score(r2ce:float, r2bm:float, r2lat:float, cviol:float) -> float:
        c,b,l,v = map(float,[r2ce,r2bm,r2lat,cviol])
        return 2*c + b + l - v/5

    # menv = Env(Import('math','log10'))

    sq   = Query({'f':Fit.id(),
                 **{x:Fit[x]() for x in resid_cols}})
    spb  = PyBlock(f_score, args = [sq[x] for x in resid_cols])
    score =                                                                     \
        Gen(name    = 'score',
            desc    = 'Computes Fit.score',
            query   = sq,
            funcs   = [spb],
            tags    = ['fit'],
            actions = [Fit(fit   = sq['f'],
                           score = spb['out'])])

    ############################################################################

    fm_env = defaultEnv + Env(Import('functionals','FromMatrix'))

    def get_lda_viol(x : str, a1: float, msb: float) -> float:
        out = abs(1 - FromMatrix(np.array(loads(x)),float(a1),float(msb)).apply(s=0,a=1))
        return out

    lvpb = PyBlock(get_lda_viol,
                   env  = fm_env,
                   args = [fbq[x] for x in 'ram'])

    lda_viol =                                                                 \
        Gen(name    = 'lda_viol',
            desc    = 'Computes Fit.lda_viol',
            query   = fbq,
            funcs   = [lvpb],
            tags    = ['fit'],
            actions = [Fit(fit=fbq['f'],lda_viol=lvpb['out'])])

    ############################################################################

    def get_h_viol(x : str, a1: float, msb: float) -> float:
        hvec = one_time_function(float(a1),float(msb))
        return abs(.3125 + hvec @ np.array(loads(x)))

    hvpb = PyBlock(get_h_viol,
                   env  = fm_env+Env(Import('functionals.scripts.fit.h_norm_const','one_time_function')),
                   args = [fbq[x] for x in 'ram'])

    h_viol =                                                                 \
        Gen(name    = 'h_viol',
            desc    = 'Computes Fit.h_viol',
            query   = fbq,
            funcs   = [hvpb],
            tags    = ['fit','parallel'],
            actions = [Fit(fit = fbq['f'], h_viol = hvpb['out'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [resid,score,lda_viol,h_viol]

    mod.add(gens)
