from typing  import Tuple as T, List as L, Optional as O
import numpy as np # type: ignore
from json import loads, dumps
from os      import environ
from os.path import join

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, COALESCE, COUNT, EQ, CONVERT, Text,
                    flatten, Env, Import, defaultEnv, Const, Constraint, JPath,
                    GROUP_CONCAT, MAX, AND, GT)

from functionals.scripts.fit.h_norm_const     import  one_time_function
from functionals.scripts.fit.opt              import  opt
from functionals.scripts.fit.db_data          import  db_data
from functionals.scripts.fit.a_ce             import a_ce
from functionals.scripts.fit.fit_vectors      import fit_vectors
from functionals.scripts.fit.mk_plt           import mk_plt

############################################################################
############################################################################
############################################################################
functionals_db = '/'+join(*__file__.split('/')[:-3],'data/functionals.json')

def fit(mod : Model) -> None:
    # Extract tables

    # Extract tables
    tabs = ['fit','calc','fitparams','functional','bulks','refs']
    Fit,Calc, Fitparams,Fx,Bulks,Refs = map(mod.get, tabs) # type: ignore

    # Extract rels
    calc__functional,fit__calc,bulks__calc,refs__bulks,\
    fit__fitparams\
     = map(mod.get_rel, [Calc.r('functional'),Fit.r('calc'),
                        Bulks.r("calc"),Refs.r('bulks'),Fit.r('fitparams')])
    ############################################################################
    ############################################################################

    refpth  = JPath("refs", [refs__bulks])
    beefpth = JPath('functional',[calc__functional,bulks__calc])
    aceq = Query(exprs=dict(b = Bulks.id(),
                            r = GROUP_CONCAT(Refs['energy'](refpth)),
                            c = GROUP_CONCAT(Refs['contribs'](refpth),delim='$'),
                            e = Bulks['eng'](),
                            x = Bulks['contribs'](),
                            n = Bulks['n_atoms'](),
                            f = MAX(Fx['data'](beefpth)),
                            m = Bulks['name']()),
                 aggcols = [Bulks.id()],
                 basis   = [Bulks])

    acepb = PyBlock(a_ce, args = [aceq[x] for x in 'rcexnfm'],)

    ace =                                                                       \
        Gen(name    = 'ace',
            desc    = 'Computes A matrix vectors to yield cohesive energy',
            query   =  aceq,
            funcs   = [acepb],
            actions = [Bulks(bulks=aceq['b'],ab_ce=acepb['out'])])

    ############################################################################
    vecout = ['ab_'+x for x in ['bm','vol']]

    aceq = Query(exprs=dict(b = Bulks.id(),
                            r = GROUP_CONCAT(Refs['energy'](refpth)),
                            c = GROUP_CONCAT(Refs['contribs'](refpth),delim='$'),
                            x = Bulks['contribs'](),
                            e = Bulks['energies'](),
                            n = Bulks['n_atoms'](),
                            f = MAX(Fx['data'](beefpth))),
                 aggcols = [Bulks.id()],
                 basis   = [Bulks])


    vecq = Query(exprs=dict(b = Bulks.id(),
                            e = Bulks['energies'](),
                            v = Bulks['volumes'](),
                            c = Bulks['contribs'](),
                            f = Fx['data'](beefpth),
                            l = Bulks['volume'](),
                            n = Bulks['name']()),
                 basis   = [Bulks],
                 opt_attr = [Bulks['expt_bm'](), Bulks['expt_vol']()])

    vecpb = PyBlock(fit_vectors,
                    args     = [vecq[x] for x in 'evcfln'],
                    outnames = vecout)

    vecs =                                                                       \
        Gen(name    = 'vecs',
            desc    = 'Computes Ax+b=y matrices/vectors',
            query   = vecq, funcs   = [vecpb],
            tags    = ['fit'],
            actions = [Bulks(bulks = vecq['b'], **{x:vecpb[x] for x in vecout})])

    ############################################################################

    fxpth = JPath('functional',[calc__functional])
    bpth  = JPath('bulks',[bulks__calc])
    fdgbs = dict(zip('abcdefgh',['ab_ce','ab_bm','ab_vol',
                                    'expt_ce','expt_bm','expt_vol','name','ce']))

    def null_to_empty(x:Expr)->Expr:
        '''Cast attribute to string - if null, return empty string.'''
        return COALESCE([CONVERT(x,Text()), Lit('')])

    fddict = {k:GROUP_CONCAT(null_to_empty(Bulks[v](bpth)),delim='$')
                for k,v in fdgbs.items()}

    fdq = Query(exprs    = dict(z = Calc.id(), **fddict),
                basis    = [Calc],
                aggcols  = [Calc.id()],
                opt_attr = [Bulks[v](bpth) for v in fdgbs.values()],
                constr   = AND([Fx['beef'](fxpth),]))
                               #Calc['done']()]))

    fdpb = PyBlock(db_data, args = [fdq[x] for x in 'abcdefgh'])

    fitdata = Gen(
        name    = 'fitdata',
        desc    = 'Assembles all BEEF data into objects for fitting',
        query   = fdq,
        funcs   = [fdpb],
        tags    = ['fit'],
        actions = [Calc(calc = fdq['z'], fitdata = fdpb['out'])]
    )

    ############################################################################

    fpth     = JPath('fitparams',[fit__fitparams])
    errs     = {k : Fitparams[v+'_scale'](fpth)
                for k,v in zip('cbl',['ce','bm','lc'])}

    pq = Query(exprs  = dict(f = Fit.id(), n=Fit['name'](),
                             p = Fit['pth'](), **errs),
               basis  = [Fit],
               constr = GT(Fit['done'](),Lit(0)))

    plt_out = ['plt','err_ce','err_bm','err_lat']
    ppb = PyBlock(mk_plt,args = [pq[x] for x in 'npcbl'],outnames=plt_out)

    plt =                                                                     \
        Gen(name    = 'plt',
            desc    = 'Plots metrics for a fit result',
            query   = pq,
            tags    = ['fit'],
            funcs   = [ppb],
            actions = [Fit(fit=pq['f'],**{x:ppb[x] for x in plt_out})])

    ############################################################################
    funmetrics = ['mse_%s'%(x) for x in ['ce','bm','lat']]

    def fefunc(pth:str)->T[float,float,float]:
        from functionals.fit.fit import FitResult
        traj = FitResult.from_pth(pth).full
        ce,bm,lat = [traj.test_err[x][-1] for x in ['ce','bm','lc']]
        return ce,bm,lat

    feq = Query(dict(f = Fit.id(), p = Fit['pth']()),constr=GT(Fit['done'](),Lit(0)))
    fepb = PyBlock(fefunc,args=[feq['p']],outnames=funmetrics)
    fiteval =                                                                     \
        Gen(name    = 'fiteval',
            desc    = 'MSE of most constrained point along fit trajectory',
            query   = feq,
            tags    = ['fit'],
            funcs   = [fepb],
            actions = [Fit(fit=feq['f'],**{x:fepb[x] for x in funmetrics})])


    ############################################################################
    # r2metrics = ['r2_'+x for x in ['ce','bm','lat']]
    # r2q = Query(exprs  = dict(f = Fit.id(), p=Fit['pth']()))
    # r2pb = PyBlock(r2_fit, args=[r2q['p']],outnames=r2metrics)
    # fitr2 = Gen(name='fit_r2',
    #            desc='populates fit.r2 values',
    #            query=r2q, tags=['fit'],funcs=[r2pb],
    #            actions=[Fit(fit=r2q['f'],**{x:r2pb[x] for x in r2metrics})])

    ############################################################################
    ############################################################################
    ############################################################################
    gens = [ace,vecs,fitdata,plt,fiteval]
    mod.add(gens)
