# External Modules
from typing import Tuple as T, List as L, Optional as O
from csv    import DictReader, reader
from math   import log10
from ast    import literal_eval
from os     import environ
from json   import loads
import numpy as np # type: ignore

# Internal Modules
from dbgen import (Model, Generator as Gen, Expr, Literal, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, AS, COALESCE, COUNT,
                    SimpleFunc,flatten, Env, Import, defaultEnv)

from functionals.fit import FromMatrix,Functional,BEEF

from functionals.scripts.load.parse_fitstep   import  parse_fitstep
from functionals.scripts.fit.fitting          import  fitting
from functionals.scripts.fit.h_norm_const     import  h_norm_const

################################################################################
################################################################################
################################################################################

functionals_db = '/Users/ksb/Documents/JSON/functionals.json'
################################################################################


################################################################################
csv_env = defaultEnv + Env(csv=['reader'])
def const() -> T[L[O[str]],L[O[str]],L[O[float]],L[O[float]],L[O[float]],L[O[str]],L[O[str]]]:

    with open('/Users/ksb/functionals/data/constraints.csv','r') as f:
        r = reader(f)
        next(r)
        name,desc,s,a,val,typ,vec = tuple(map(list,zip(*r)))
    s_,a_,val_  = [[float(str(x)) if x else None for x in xs] for xs in [s,a,val]]
    n_,d_,t_,v_ = [[str(x) if x else None for x in xs] for xs in [name,desc,typ,vec]]
    return n_,d_,s_,a_,val_,t_,v_

def nlconst() -> T[L[O[str]],L[O[str]],L[O[str]],L[O[str]],L[O[str]],L[O[float]],L[O[float]]]:

    with open('/Users/ksb/functionals/data/nlconstraints.csv','r') as fi:
        r = reader(fi)
        next(r)
        name,desc,f,j,h,lb,ub = tuple(map(list,zip(*r)))
    l_,u_  = [[float(str(x)) if x else None for x in xs] for xs in [lb,ub]]
    n_,d_,f_,j_,h_ = [[str(x) if x else None for x in xs] for xs in [name,desc,f,j,h]]
    return n_, d_, f_, j_, h_, l_, u_

################################################################################
fitcols = ['name','consts','nlconsts','dataconst','basis','initfit','bound',
           'maxiter','constden','bm_weight','lat_weight']

pf_env = Env(os = ['environ'], csv = ['DictReader']) + defaultEnv

def parse_fits() -> T[L[str]]:
    fitcols = ['name','consts','nlconsts','dataconst','basis','initfit','bound',
               'maxiter','constden','bm_weight','lat_weight']
    vals = {c:[] for c in fitcols} # type: dict
    with open(environ['FITPATH'],'r') as f:
        reader = DictReader(f)
        for row in reader:
            for k,v in vals.items(): v.append(row[k])
    return tuple(map(list,vals.values())) # type: ignore

################################################################################
def sqldict(d : dict) -> Expr:
    '''Turn a Python dict into a SQL (serialized) dict'''
    args = [[',', ' "%s" : ['%k, GC(v), ']'] for k,v in d.items()]
    return C('{',*flatten(args)[1:],'}')

################################################################################
################################################################################
################################################################################

def fit(mod:Model) -> None:
    # Extract tables
    tabs = ['fit','const','fit_step','fit_const','fit_data','dft_data',
            'nonlin_const','fit_nonlin_const','globals','expt','species']

    Fit, Const, Fit_step, Fit_const, Fit_data, D, Nl_const, Fit_nl_const,\
    Globals,Expt,Species = map(mod.get, tabs) # type: ignore

    ########################################################################

    constcols = ['const_name','description','s','alpha','val','kind','vec']
    pop_constraint =                                                            \
        Gen(name    = 'pop_constraint',
            desc    = 'populate constraints',
            targets = [Const,*Const.select(constcols[1:])],
            func    = SimpleFunc((const,csv_env), outputs = constcols),
            tags    = 'fit')

    ########################################################################
    nlconstcols = ['nlconst_name','description','f','df','hess','lb','ub']
    pop_nlconstraint =                                                          \
        Gen(name    = 'pop_nlconstraint',
            desc    = 'populate nonlinear constraints',
            targets = [Nl_const,*Nl_const.select(nlconstcols[1:])],
            func    = SimpleFunc((nlconst,csv_env), outputs = nlconstcols),
            tags    = 'fit')

    ########################################################################
    all_fit = dict(name='all',consts = '{}',nlconsts = '{}', dataconst='.',
                   basis = 8, initfit = 1, bound = 1, maxiter=1000,
                   constden = 3, lat_weight = 1, bm_weight = 1)
    pop_all_fit =                                                           \
        Gen(name    = 'pop_all_fit',
            desc    = 'Fills a special row in Fit table that includes all data and constraints',
            targets = [Fit, *Fit.select(fitcols[1:])],
            consts  = all_fit,
            tags    = 'fit')
    ########################################################################
    pop_fit =                                                               \
        Gen(name    = 'pop_fit',
            desc    = 'Looks in FITPATH for fitting specifications',
            targets = [Fit, *Fit.select(fitcols[1:])],
            func    = SimpleFunc((parse_fits,pf_env),outputs=fitcols),
            tags    = ['fit','io'])
    ########################################################################

    pop_fitdata =                                                           \
        Gen(name    = 'pop_fitdata',
            desc    = 'Connect fit to DFT data based on Regex match',
            targets = Fit_data,
            basis   = ['fit','dft_data'],
            constr  = REGEXP(Species.nickname,BINARY(Fit.dataconst)),
            tags    = 'fit')

    ########################################################################

    def parse_dict(dstr:str,cname:str)->T[str,float]:

        if not dstr:
            return [],[] # type: ignore
        else:
            d = literal_eval(dstr)
            return (cname, d.get(cname,1.0))

    pfc = defaultEnv + Env(ast=['literal_eval'])

    pop_fitconst =                                                          \
        Gen(name    = 'pop_fitconst',
            desc    = 'Parses weight dictionary',
            targets = [Fit_const, Fit_const.const_weight],
            input   = [Fit.consts, Const.const_name |AS| 'cname'],
            basis   = ['Fit','Const'],
            func    = SimpleFunc((parse_dict, pfc),
                                    ['consts','cname'],
                                    ['const_name','const_weight']),
            tags    =  'fit')

    pop_fitnlconst =                                                            \
        Gen(name    = 'pop_fitnlconst',
            desc    = 'Parses weight dictionary' ,
            targets = [Fit_nl_const, Fit_nl_const.nl_const_weight],
            input   = [Fit.nlconsts, Nl_const.nlconst_name |AS| 'cname'],
            basis   = ['fit',Nl_const._name],
            func    =  SimpleFunc((parse_dict, pfc),
                                    ['nlconsts','cname'],
                                    ['nlconst_name','nl_const_weight']),
            tags    = 'fit')

    ########################################################################

    pfs_env = defaultEnv + Env(re=['findall'])

    pop_fitstep =                                                               \
        Gen(name    = 'pop_fitstep',
            desc    = 'Analyzes fitting result log, populates Fitstep',
            targets = [Fit_step, Fit_step.cost, Fit_step.c_viol],
            input   =  Fit.log,
            func    = SimpleFunc((parse_fitstep, pfs_env),
                                  outputs=['niter','cost','c_viol']),
            tags    = ['fit', 'parallel'])

    ########################################################################

    data_dict = dict(
        atomic_contribs      = C('"',D.atomic_contribs,'"'),
        atomic_energies      = C('"',D.atomic_energies,'"'),
        bulk_ratio           = D.bulk_ratio,
        bulk_energy          = D.bulk_energy,
        bulk_contribs        = C('"',D.bulk_contribs,'"'),
        coefs                = C('"',D.coefs,'"'),
        composition          = C('"',D.composition,'"'),
        expt_cohesive_energy = D.expt_cohesive_energy,
        expt_volume          = D.expt_volume,
        expt_bm              = D.expt_bm,
        energy_vector        = D.energy_vector,
        volume_vector        = D.volume_vector,
        contrib_vector       = D.contrib_vector)

    raw_data =                                                              \
        Gen(name    = 'fit_rawdata',
            desc    = 'Serialize info stored in dft_data mapping table',
            targets = Fit.raw_data,
            input   = (sqldict(data_dict) |AS| 'raw_data'),
            links   = Fit_data.Fit,
            agg     = 'fit',
            tags    = 'fit')

    ########################################################################

    constdict = dict(name   = C('"',Const.const_name,'"'),
                     s      = COALESCE(Const.s,'None'),
                     alpha  = COALESCE(Const.alpha,'None'),
                     val    = COALESCE(Const.val,'None'),
                     kind   = C('"',Const.kind,'"'),
                     vec    = C('"',Const.vec,'"'),
                     weight = Fit_const.const_weight)

    raw_const =                                                                 \
        Gen(name    = 'raw_const',
            desc    =  'Serialize info stored in constraint mapping table',
            targets = Fit.raw_const,
            input   = sqldict(constdict) |AS| 'raw_const',
            links   = Fit_const.Fit,
            agg     = 'fit',
            option  = [Const.s,Const.alpha],
            tags    = 'fit')

    ########################################################################

    nlcdict = dict(name   = C('"',Nl_const.nlconst_name,'"'),
                   f      = C('"',Nl_const.f,'"'),
                   df     = C('"',Nl_const.df,'"'),
                   hess   = C('"',Nl_const.hess,'"'),
                   ub     = Nl_const.ub,
                   lb     = Nl_const.lb,
                   weight = Fit_nl_const.nl_const_weight)

    raw_nlconst =                                                           \
        Gen(name    = 'raw_nlconst',
            desc    = 'Serialize info stored in nl_const mapping table',
            targets = Fit.raw_nlconst,
            input   = sqldict(nlcdict) |AS| 'raw_nlconst',
            links   = Fit_nl_const.Fit,
            agg     = 'fit',
            tags    ='fit')

    ########################################################################

    fit_out = ['log','timestamp','runtime','result','err']

    fit_in  = ['raw_data','raw_const','raw_nlconst','maxiter','basis','bound',
                'initfit','constden','bm_weight','lat_weight']

    fit_env = Env([Import('functionals.fit',objs=['Fit'])]) + defaultEnv

    make_fit =                                                              \
        Gen(name    = 'make_fit',
            desc    = "Scipy minimize: fill out 'Fit.err' iff it fails",
            targets = Fit.select(fit_out),
            input   = Fit.select(fit_in),
            func    = SimpleFunc((fitting, fit_env),
                                   inputs  = fit_in,
                                   outputs = fit_out),
            option  = [Fit.raw_data, Fit.raw_const,Fit.raw_nlconst],
            tags    = ['fit','parallel'])

    ########################################################################

    countdict = {Fit.n_const   : Fit_const,
                 Fit.n_nlconst : Fit_nl_const,
                 Fit.n_data    : Fit_data}

    nconst,nnlconst,ndata =  [                                                 \
        Gen(name    = x.name,
            desc    = 'Aggregate and count Fit-related mapping tables',
            targets = x ,
            input   = COUNT(Literal(1)),
            basis   = y._name,
            agg     = 'fit',
            tags    = 'fit')
            for x,y in countdict.items()]

    ########################################################################
    gcols  = ['all_data','all_constraints','all_nlconstraints']
    fcols  = ['raw_data','raw_const','raw_nlconst']
    zipped = zip(fcols,gcols)

    globals =                                                               \
        Gen(name    = 'globals',
            desc    = 'Copies data from "all" row in Fit to Globals table',
            targets = [Globals,*Globals.select(gcols)],
            input   = [Fit._get(x)|AS|y for x,y in zipped],
            constr  = Fit.name == 'all',
            tags    = 'fit')

    ########################################################################
    resid_cols = ['r2_ce','r2_bm','r2_lat','c_viol']

    fm_env = defaultEnv + Env(functionals=['FromMatrix'])

    def f(dat:str,con:str,nlcon:str,x:str,n:int)->T[float,float,float,float]:
        return FromMatrix(np.array(loads(x)).reshape((n,n))).costs(dat,con,nlcon)

    resid =                                                                 \
        Gen(name    = 'resid',
            desc    = 'Computes metrics for a fit result',
            targets = Fit.select(resid_cols),
            input   = Globals.select(gcols)
                                  +[Fit.result,Fit.basis],
            basis   = ['Fit','Globals'],
            tags    = 'fit',
            func    = SimpleFunc((f,fm_env),
                                 inputs= gcols+['result','basis'],
                                 outputs=resid_cols))
    ########################################################################

    def bdist(x:str,n:int)->float:
        fx = FromMatrix(np.array(loads(x)).reshape((n,n)))
        r  = np.logspace(-2,2,5).tolist()
        return Functional.diff(fx,BEEF,r,r)

    beef_env = Env(functionals=['BEEF','Functional'])

    beefdist =                                                              \
        Gen(name = 'beefdist',
            desc = 'Computes Fit.beefdist',
            targets = Fit.beefdist,
            input = [Fit.result,Fit.basis],
            func = (bdist, fm_env + beef_env),
            tags = 'fit')

    ########################################################################

    def f_score(r2ce:float, r2bm:float, r2lat:float, cviol:float) -> float:
        return float(2*r2ce + r2bm + r2lat) - log10(cviol)/5

    menv = Env(math=['log10'])

    score =                                                                 \
        Gen(name    = 'score',
            desc    = 'Computes Fit.score',
            targets = Fit.score,
            input   = Fit.select(resid_cols),
            func    = SimpleFunc((f_score, menv),
                                  inputs = resid_cols, outputs = ['score']),
            tags    = 'fit')

    ########################################################################

    def get_lda_viol(n : int, x : str) -> float:
        M = np.array(loads(x)).reshape((n,n))
        return abs(1 - FromMatrix(M).apply(s=0,a=1))

    lda_viol =                                                                 \
        Gen(name    = 'lda_viol',
            desc    = 'Computes Fit.lda_viol',
            targets = Fit.lda_viol ,
            input   = [Fit.result, Fit.basis],
            func    = (SimpleFunc((get_lda_viol,fm_env),
                                                ['basis','result'],
                                                ['lda_viol'])),
            tags    = 'fit')
    gens = [pop_constraint,pop_nlconstraint,pop_all_fit,pop_fit,pop_fitdata,
            pop_fitconst,pop_fitnlconst,pop_fitstep,raw_data,raw_const,
            raw_nlconst,nconst,nnlconst,ndata,globals,resid,beefdist,score,
            lda_viol]
    mod.add(gens)
