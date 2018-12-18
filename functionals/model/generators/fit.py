# External Modules
from typing import Tuple as T, List as L, Optional as O
from csv    import DictReader, reader
from math   import log10
from ast    import literal_eval
from os     import environ
from json   import loads
import numpy as np # type: ignore

# Internal Modules
from dbgen2 import (Model, Gen, Ref, Query, PyBlock, Expr, Literal, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, COALESCE, COUNT, EQ,
                    flatten, Env, Import, defaultEnv, Const as Constant)

from functionals.fit import FromMatrix,Functional,BEEF

from functionals.scripts.load.parse_fitstep   import  parse_fitstep
from functionals.scripts.fit.fitting          import  fitting
from functionals.scripts.fit.h_norm_const     import  h_norm_const

################################################################################
################################################################################
################################################################################

functionals_db = '/Users/ksb/Documents/JSON/functionals.json'

################################################################################
csv_env = defaultEnv + Env(Import('csv','reader'))
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

pf_env = Env(Import('os','environ'), Import('csv','DictReader')) + defaultEnv

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

def fit(mod : Model) -> None:
    # Extract tables
    tabs = ['fit','const','fit_step','fit_const','fit_data','dft_data',
            'nonlin_const','fit_nonlin_const','globals','expt','species']

    Fit, Const, Fit_step, Fit_const, Fit_data, D, Nl_const, Fit_nl_const,\
    Globals,Expt,Species = map(mod.get, tabs) # type: ignore


    fbq = Query({'r':Fit['result'],'b':Fit['basis']})

    ########################################################################

    constcols = ['const_name','description','s','alpha','val','kind','vec']
    pcpb = PyBlock(const, env = csv_env, outnames = constcols)
    pop_constraint =                                                            \
        Gen(name    = 'pop_constraint',
            desc    = 'populate constraints',
            actions = [Const(insert = True, **{x:pcpb[x] for x in constcols})],
            funcs   = [pcpb],
            tags    = ['fit'])

    ########################################################################
    nlconstcols = ['nlconst_name','description','f','df','hess','lb','ub']
    ncpb = PyBlock(nlconst,env = csv_env, outnames = nlconstcols)
    pop_nlconstraint =                                                          \
        Gen(name    = 'pop_nlconstraint',
            desc    = 'populate nonlinear constraints',
            actions = [Nl_const(insert = True, **{x:ncpb[x] for x in nlconstcols})],
            funcs   = [ncpb],
            tags    = ['fit'])

    ########################################################################
    all_fit = dict(name='all',consts = '{}',nlconsts = '{}', dataconst='.',
                   basis = 8, initfit = 1, bound = 1, maxiter=1000,
                   constden = 3, lat_weight = 1, bm_weight = 1)
    afpb = PyBlock(lambda x:x, args = [Constant(tuple(all_fit.values()))],
                    outnames = list(all_fit.keys()))
    pop_all_fit =                                                               \
        Gen(name    = 'pop_all_fit',
            desc    = 'Fills a special row in Fit table that includes all data and constraints',
            actions = [Fit(insert = True, **{x:afpb[x] for x in all_fit.keys()})],
            funcs   = [afpb],
            tags    = ['fit'])
    ########################################################################
    pfpb = PyBlock(parse_fits,env = pf_env,outnames=fitcols)
    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'Looks in FITPATH for fitting specifications',
            actions = [Fit(insert = True, **{x:pfpb[x] for x in fitcols})],
            funcs   = [pfpb],
            tags    = ['fit','io'])
    ########################################################################
    pfdq = Query(exprs   = {},
                 basis   = ['fit','dft_data'],
                 constr  = REGEXP(Species['nickname'],BINARY(Fit['dataconst'])))
    pop_fitdata =                                                               \
        Gen(name    = 'pop_fitdata',
            desc    = 'Connect fit to DFT data based on Regex match',
            actions = [Fit_data(insert   = True,
                                fit      = Ref('fit'),
                                dft_data = Ref('dft_data'))],
            query   = pfdq,
            tags    = ['fit'])

    ########################################################################

    def parse_dict(dstr:str,cname:str)->T[str,float]:

        if not dstr:
            return [],[] # type: ignore
        else:
            d = literal_eval(dstr)
            return (cname, d.get(cname,1.0))

    pfc = defaultEnv + Env(Import('ast','literal_eval'))

    pfcq = Query({'c':Fit['consts'],'cn':Const['const_name']},
                  basis   = ['Fit','Const'])

    pfcpb = PyBlock(parse_dict, env = pfc,
                    args = [pfcq['c'],pfcq['cn']],
                    outnames = ['cn','cw'])
    pop_fitconst =                                                          \
        Gen(name    = 'pop_fitconst',
            desc    = 'Parses weight dictionary',
            actions = [Fit_const(insert = True,
                                 fit    = Ref('fit'),
                                 const  = Const(const_name=pfcpb['cn']),
                                 const_weight = pfcpb['cw'])],
            query   = pfcq,
            funcs   = [pfcpb],
            tags    = ['fit'])

    pfncq = Query({'c':Fit['nlconsts'],'n':Nl_const['nlconst_name']},
                  basis   = ['Fit','nonlin_const'])

    pfncpb = PyBlock(parse_dict, env = pfc,
                    args = [pfncq[x] for x in 'cn'],
                    outnames = ['n','w'])
    pop_fitnlconst =                                                            \
        Gen(name    = 'pop_fitnlconst',
            desc    = 'Parses weight dictionary' ,
            actions = [Fit_nl_const(insert = True,
                                    fit    = Ref('fit'),
                                    nl_const_weight = pfncpb['w'],
                                    nonlin_const = Nl_const(nlconst_name=pfncpb['n']))],
            query   = pfncq,
            funcs   =  [pfncpb],
            tags    = ['fit'])

    ########################################################################

    pfs_env = defaultEnv + Env(Import('re','findall'))
    pfsq    = Query({'log':Fit['log']})
    pfsns   = ['niter','cost','c_viol']
    pfspb   = PyBlock(parse_fitstep, env = pfs_env,
                      args=[pfsq['log']], outnames = pfsns)
    pop_fitstep =                                                      \
        Gen(name    = 'pop_fitstep',
            desc    = 'Analyzes fitting result log, populates Fitstep',
            actions = [Fit_step(insert = True,
                                fit    = Ref('fit'),
                                **{x:pfspb[x] for x in pfsns})],
            query   =  pfsq,
            funcs   = [pfspb],
            tags    = ['fit', 'parallel'])

    ########################################################################

    data_dict = dict(
        atomic_contribs      = C('"',D['atomic_contribs'],'"'),
        atomic_energies      = C('"',D['atomic_energies'],'"'),
        bulk_ratio           = D['bulk_ratio'],
        bulk_energy          = D['bulk_energy'],
        bulk_contribs        = C('"',D['bulk_contribs'],'"'),
        coefs                = C('"',D['coefs'],'"'),
        composition          = C('"',D['composition'],'"'),
        expt_cohesive_energy = D['expt_cohesive_energy'],
        expt_volume          = D['expt_volume'],
        expt_bm              = D['expt_bm'],
        energy_vector        = D['energy_vector'],
        volume_vector        = D['volume_vector'],
        contrib_vector       = D['contrib_vector'])

    rdq = Query({'rd':sqldict(data_dict)},
                links = [Fit_data.r('Fit')],
                basis = ['fit'],
                aggcols = [Fit['name']])
    raw_data =                                                              \
        Gen(name    = 'fit_rawdata',
            desc    = 'Serialize info stored in dft_data mapping table',
            actions = [Fit(fit=Ref('fit'),raw_data=rdq['rd'])],
            query   = rdq,
            tags    = ['fit'])

    ########################################################################

    constdict = dict(name   = C('"',Const['const_name'],'"'),
                     s      = COALESCE(Const['s'],'None'),
                     alpha  = COALESCE(Const['alpha'],'None'),
                     val    = COALESCE(Const['val'],'None'),
                     kind   = C('"',Const['kind'],'"'),
                     vec    = C('"',Const['vec'],'"'),
                     weight = Fit_const['const_weight'])


    rcq = Query({'rc':sqldict(constdict)},
                links    = [Fit_const.r('Fit')],
                opt_attr = [Const['s'],Const['alpha']],
                aggcols  = [Fit['name']])
    raw_const =                                                                 \
        Gen(name    = 'raw_const',
            desc    =  'Serialize info stored in constraint mapping table',
            actions = [Fit(fit=Ref('fit'),raw_const=rcq['rc'])],
            query   =  rcq,
            tags    = ['fit'])

    ########################################################################

    nlcdict = dict(name   = C('"',Nl_const['nlconst_name'],'"'),
                   f      = C('"',Nl_const['f'],'"'),
                   df     = C('"',Nl_const['df'],'"'),
                   hess   = C('"',Nl_const['hess'],'"'),
                   ub     = Nl_const['ub'],
                   lb     = Nl_const['lb'],
                   weight = Fit_nl_const['nl_const_weight'])

    rncq = Query({'rnc':sqldict(nlcdict)},
                links    = [Fit_nl_const.r('Fit')],
                aggcols  = [Fit['name']])
    raw_nlconst =                                                           \
        Gen(name    = 'raw_nlconst',
            desc    = 'Serialize info stored in nl_const mapping table',
            actions = [Fit(fit=Ref('fit'),raw_nlconst=rncq['rnc'])],
            query   = rncq,
            tags    =['fit'])

    ########################################################################

    fit_out = ['log','timestamp','runtime','result','err']

    fit_in  = ['raw_data','raw_const','raw_nlconst','maxiter','basis','bound',
                'initfit','constden','bm_weight','lat_weight']

    fit_env = Env(Import('functionals.fit','Fit')) + defaultEnv

    mfq = Query({x:Fit[x] for x in fit_in},
                opt_attr  = [Fit['raw_data'],
                             Fit['raw_const'],
                             Fit['raw_nlconst']])

    mfpb = PyBlock(fitting, env = fit_env,
                   args  = [mfq[x] for x in fit_in],
                   outnames = fit_out)
    make_fit =                                                              \
        Gen(name    = 'make_fit',
            desc    = "Scipy minimize: fill out 'Fit['err' iff it fails",
            actions = [Fit(fit=Ref('fit'),**{x:mfpb[x] for x in fit_out})],
            query   = mfq,
            funcs   = [mfpb],
            tags    = ['fit','parallel'])

    ########################################################################

    countdict = {'n_const'   : Fit_const,
                 'n_nlconst' : Fit_nl_const,
                 'n_data'    : Fit_data}

    threegens = []

    for x,y in countdict.items():
        q = Query({'n':COUNT(Literal(1))},
                  aggcols = [Fit.id()],
                  basis = [y.name]
                  )
        threegens.append(
            Gen(name    = 'pop_'+x,
                desc    = 'Aggregate and count Fit-related mapping tables',
                actions = [Fit(fit=Ref('fit'),**{x:q['n']})],
                query   = q,
                tags    = ['fit']))

    nconst,nnlconst,ndata = threegens
    ########################################################################
    gcols  = ['all_data','all_constraints','all_nlconstraints']
    fcols  = ['raw_data','raw_const','raw_nlconst']
    zipped = zip(fcols,gcols)

    gq = Query({y:Fit[x] for x,y in zipped},
                constr  = EQ(Fit['name'], 'all'))


    globals =                                                               \
        Gen(name    = 'globals',
            desc    = 'Copies data from "all" row in Fit to Globals table',
            actions = [Globals(insert = True, **{g:gq[g] for g in gcols})],
            query   = gq,
            tags    = ['fit'])

    ########################################################################
    resid_cols = ['r2_ce','r2_bm','r2_lat','c_viol']

    fm_env = defaultEnv + Env(Import('functionals','FromMatrix'))

    def f(dat:str,con:str,nlcon:str,x:str,n:int)->T[float,float,float,float]:
        return FromMatrix(np.array(loads(x)).reshape((n,n))).costs(dat,con,nlcon)

    rq = Query({'res':Fit['result'],'b':Fit['basis'],
                **{g:Globals[g] for g in gcols}},
                basis   = ['Fit','Globals'])
    rpb = PyBlock(f, env = fm_env,
                  args= [rq[x] for x in gcols+['res','b']],
                  outnames = resid_cols)
    resid =                                                                     \
        Gen(name    = 'resid',
            desc    = 'Computes metrics for a fit result',
            actions = [Fit(fit=Ref('fit'),**{x:rpb[x] for x in resid_cols})],
            query   = rq,
            tags    = ['fit'],
            funcs   = [rpb])
    ########################################################################

    def bdist(x:str,n:int)->float:
        fx = FromMatrix(np.array(loads(x)).reshape((n,n)))
        r  = np.logspace(-2,2,5).tolist()
        return Functional.diff(fx,BEEF,r,r)

    beef_env = Env(Import('functionals','BEEF','Functional'))

    bdpb = PyBlock(bdist, env = fm_env + beef_env,
                    args=[fbq['r'],fbq['b']], outnames = ['d'])
    beefdist =                                                                  \
        Gen(name    = 'beefdist',
            desc    = 'Computes Fit.beefdist',
            actions = [Fit(fit=Ref('fit'),beefdist=bdpb['d'])],
            query   = fbq,
            funcs   = [bdpb],
            tags    = ['fit'])

    ########################################################################

    def f_score(r2ce:float, r2bm:float, r2lat:float, cviol:float) -> float:
        return float(2*r2ce + r2bm + r2lat) - log10(cviol)/5

    menv = Env(Import('math','log10'))
    sq   = Query({x:Fit[x] for x in resid_cols})
    spb  = PyBlock(f_score, env = menv,
                  args = [sq[x] for x in resid_cols], outnames = ['s'])
    score =                                                                     \
        Gen(name    = 'score',
            desc    = 'Computes Fit.score',
            actions = [Fit(fit=Ref('fit'),score=spb['s'])],
            query   = sq,
            funcs   = [spb],
            tags    = ['fit'])

    ########################################################################

    def get_lda_viol(n : int, x : str) -> float:
        M = np.array(loads(x)).reshape((n,n))
        return abs(1 - FromMatrix(M).apply(s=0,a=1))

    lvpb = PyBlock(get_lda_viol, env = fm_env,
                    args=[fbq['r'],fbq['b']],
                    outnames = ['lv'])
    lda_viol =                                                                 \
        Gen(name    = 'lda_viol',
            desc    = 'Computes Fit.lda_viol',
            actions = [Fit(fit=Ref('fit'),lda_viol=lvpb['lv'])] ,
            query   = fbq,
            funcs   = [lvpb],
            tags    = ['fit'])


    ########################################################################
    ########################################################################
    ########################################################################

    gens = [pop_constraint,pop_nlconstraint,pop_all_fit,pop_fit,pop_fitdata,
            pop_fitconst,pop_fitnlconst,pop_fitstep,raw_data,raw_const,
            raw_nlconst,nconst,nnlconst,ndata,globals,resid,beefdist,score,
            lda_viol, make_fit]

    mod.add(gens)
