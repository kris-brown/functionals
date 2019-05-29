from typing  import Tuple as T, List as L, Optional as O
import numpy as np # type: ignore
from json import loads, dumps
from os      import environ
from os.path import join

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, COALESCE, COUNT, EQ, CONVERT, Text,
                    flatten, Env, Import, defaultEnv, Const, Constraint, JPath,
                    GROUP_CONCAT)

from functionals.scripts.fit.h_norm_const     import  one_time_function
from functionals.scripts.fit.opt              import  opt
from functionals.scripts.fit.db_data          import  db_data

############################################################################
############################################################################
############################################################################
functionals_db = join(environ['FUNCTIONALS_ROOT'],'data/functionals.json')

def fit(mod : Model) -> None:
    # Extract tables

    # Extract tables
    tabs = ['fit','job','calc','fitparams','functional','beef','bulks']
    Fit, Job,Calc, Fitparams,Fx,BEEF,Bulks = map(mod.get, tabs) # type: ignore

    # Extract rels
    beef__functional,calc__functional,fit__calc,bulks__job,job__calc = map(mod.get_rel,
        [BEEF.r('functional'),Calc.r('functional'),Fit.r('calc'),Bulks.r("job"),Job.r('calc')])
    ############################################################################
    ############################################################################
    ############################################################################
    fxpth = JPath('functional',[calc__functional])
    bpth = JPath('bulks',[bulks__job,job__calc])
    fdgbs = dict(zip('abcdefghijk',['a_ce','a_bm','a_l','b_ce','b_bm','b_l',
                                    'expt_ce','expt_bm','expt_vol','name','ce']))
    fddict = {k:GROUP_CONCAT(COALESCE(CONVERT(Bulks[v](bpth),Text()),Lit('')),delim='$')
                for k,v in fdgbs.items()}
    fdq = Query(exprs = dict(z=Calc.id(),**fddict),
                basis = [Calc],
                aggcols= [Calc.id()],
                opt_attr = [Bulks[v](bpth) for v in fdgbs.values()],
                constr = Fx['beef'](fxpth))
    fdpb = PyBlock(db_data,args=[fdq[x] for x in 'abcdefghijk'])
    fitdata = Gen(
        name = 'fitdata',
        desc = 'Assembles all BEEF data into objects for fitting',
        query = fdq,
        funcs = [fdpb],
        actions = [Calc(calc=fdq['z'],fitdata=fdpb['out'])]
    )
    ############################################################################
    resid_cols = ['result','decaycosts','c_viol','mse_ce','mse_bm','mse_lat']


    rq = Query(dict(f = Fit.id(),
                    **{'t%d'%i : Fit['traj%d'%i]() for i in range(5)}))


    def rf(*ts:str) -> T[str,str,float,float,float,float]:
        res,dc = [],[]
        for i,t in enumerate(map(loads,ts)):
            o = t['opt']
            res.append(t['x'][o])
            dc.append(t['cost'][o])
            if i == 2:
                cviol = t['cviol_all'][o]
                msece = t['ce'][o]
                msebm = t['bm'][o]
                mselc = t['lc'][o]
        return dumps(res),dumps(dc),cviol,msece,msebm,mselc

    rpb = PyBlock(rf,
                  # env  = f_env,
                   args = [rq['t%d'%i] for i in range(5)],
                   outnames = resid_cols)

    resid =                                                                     \
        Gen(name    = 'resid',
            desc    = 'Computes metrics for a fit result',
            query   = rq,
            tags    = ['fit'],
            funcs   = [rpb],
            actions = [Fit(fit=rq['f'],
                           **{x:rpb[x] for x in resid_cols})])


    ############################################################################
    ############################################################################
    ############################################################################
    gens =  [resid,fitdata]
    mod.add(gens)
