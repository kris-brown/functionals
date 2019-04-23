from typing  import Tuple as T, List as L, Optional as O
import numpy as np # type: ignore
from json import loads, dumps
from os      import environ
from os.path import join

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, COALESCE, COUNT, EQ,
                    flatten, Env, Import, defaultEnv, Const, Constraint, JPath)

from functionals.fit import Fit as FitObj , FromMatrix #Functional,BEEF
from functionals.scripts.fit.h_norm_const     import  one_time_function
from functionals.scripts.fit.opt              import  opt

############################################################################
############################################################################
############################################################################
functionals_db = join(environ['FUNCTIONALS_ROOT'],'data/functionals.json')

def fit(mod : Model) -> None:
    # Extract tables

    # Extract tables
    tabs = ['fit','calc','fitparams','beef']
    Fit, Calc, Fitparams,BEEF = map(mod.get, tabs) # type: ignore

    # Extract rels
    beef__functional,calc__functional,fit__calc = map(mod.get_rel,
        [BEEF.r('functional'),Calc.r('functional'),Fit.r('calc')])
    ############################################################################
    ############################################################################
    ############################################################################
    resid_cols = ['decaycosts','mse_ce','mse_bm','mse_lat','r2_ce','r2_bm','r2_lat']
    f_env      = defaultEnv + Env(Import('functionals.fit.fit',Fit='FitObj'))

    # def f_resid(db:str,x:str,c:int)->T[str,float,float,float,float,float,float,float]:
    #     return db=db,x=x,calc=c)

    fc = JPath("calc", [fit__calc])

    rq = Query(exprs = dict(f = Fit.id(),
                            p = Fit['pth'](),
                            x = Fit['result']()),
                #constr = Fit['name']() |EQ| Lit('test'),
               basis = ['fit'])

    def rf(p : str, x : str) -> T[str,float,float,float,float,float,float,]:
        return FitObj.from_json(p).allcosts(x)

    rpb = PyBlock(rf,
                   env  = f_env,
                   args = [rq['p'],rq['x']],
                   outnames = resid_cols)

    resid =                                                                     \
        Gen(name    = 'resid',
            desc    = 'Computes metrics for a fit result',
            query   = rq,
            tags    = ['fit'],
            funcs   = [rpb],
            actions = [Fit(fit=rq['f'],
                           **{x:rpb[x] for x in resid_cols})])

    ########################################################################
    oq = Query(exprs = dict(f = Fit.id(),
                            p = Fit['pth']()))
    opb = PyBlock(opt,args=[oq['p']],outnames=['opt','result'])
    optg =                                                                     \
        Gen(name    = 'opt',
            query   = oq,
            funcs   = [opb],
            tags    = ['fit'],
            actions = [Fit(fit   = oq['f'],
                           opt = opb['opt'],result=opb['result'])])

    ########################################################################
    # Common queries

    def cf(p : str, x : str) -> str:
        x_ = np.array(loads(x)[2])
        return dumps(FitObj.from_json(p).allviols(x_,2))

    cvpb = PyBlock(cf, env  = f_env, args=[rq[x] for x in 'px'])

    cviol =                                                                     \
        Gen(name    = 'cviol',
            query   = rq,
            funcs   = [cvpb],
            tags    = ['fit'],
            actions = [Fit(fit = rq['f'], c_viol = cvpb['out'])])

    ########################################################################
    def f_score(r2ce:float, r2bm:float, r2lat:float, cviol:float) -> float:
        c,b,l,v = map(float,[r2ce,r2bm,r2lat,cviol])
        return 2*c + b + l - v/5

    sq   = Query({'f':Fit.id(),
                 **{x:Fit[x]() for x in resid_cols[-4:]}})
    spb  = PyBlock(f_score, args = [sq[x] for x in resid_cols[-4:]])
    score =                                                                     \
        Gen(name    = 'score',
            desc    = 'Computes Fit.score',
            query   = sq,
            funcs   = [spb],
            tags    = ['fit'],
            actions = [Fit(fit   = sq['f'],
                           score = spb['out'])])
    ############################################################################

    ############################################################################
    ############################################################################
    ############################################################################
    gens = [optg,score,resid,cviol]
    mod.add(gens)
