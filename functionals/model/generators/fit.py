# External Modules
from typing import Tuple as T, List as L, Optional as O
from csv    import DictReader, reader
from math   import log10
from ast    import literal_eval
from os     import environ
from json   import loads
import numpy as np # type: ignore

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal, REGEXP, BINARY,
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
    tabs = ['fit','const','fit_step','fit_const','fit_data',
            'nonlin_const','fit_nonlin_const','globals','expt','species']

    Fit, Const, Fit_step, Fit_const, Fit_data, Nl_const, Fit_nl_const,\
    Globals,Expt,Species = map(mod.get, tabs) # type: ignore


    fbq = Query({'f' : Fit.id, 'r' : Fit['result'], 'b' : Fit['basis']})

    ########################################################################

    constcols = ['const_name','description','s','alpha','val','kind','vec']
    pcpb = PyBlock(const, env = csv_env, outnames = constcols)
    pop_constraint =                                                            \
        Gen(name    = 'pop_constraint',
            desc    = 'populate constraints',
            funcs   = [pcpb],
            tags    = ['fit'],
            actions = [Const(insert = True, **{x:pcpb[x] for x in constcols})])

    ########################################################################
    nlconstcols = ['nlconst_name','description','f','df','hess','lb','ub']

    ncpb = PyBlock(nlconst,
                   env      = csv_env,
                   outnames = nlconstcols)
    pop_nlconstraint =                                                          \
        Gen(name    = 'pop_nlconstraint',
            desc    = 'populate nonlinear constraints',
            funcs   = [ncpb],
            tags    = ['fit'],
            actions = [Nl_const(insert = True,
                                **{x:ncpb[x] for x in nlconstcols})])

    ########################################################################
    all_fit = dict(name='all',consts = '{}',nlconsts = '{}', dataconst='.',
                   basis = 8, initfit = 1, bound = 1, maxiter=1000,
                   constden = 3, lat_weight = 1, bm_weight = 1)

    afpb = PyBlock(lambda x : x,
                   args     = [Constant(tuple(all_fit.values()))],
                   outnames = list(all_fit.keys()))

    pop_all_fit =                                                               \
        Gen(name    = 'pop_all_fit',
            desc    = 'Fills a special row in Fit table that includes all data and constraints',
            funcs   = [afpb],
            tags    = ['fit'],
            actions = [Fit(insert = True, **{x:afpb[x] for x in all_fit.keys()})])

    ########################################################################
    pfpb = PyBlock(parse_fits,
                   env      = pf_env,
                   outnames = fitcols)
    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'Looks in FITPATH for fitting specifications',
            funcs   = [pfpb],
            tags    = ['fit', 'io'],
            actions = [Fit(insert = True,
                           **{x:pfpb[x] for x in fitcols})])
    ########################################################################
    pfdq = Query(exprs   = {'f' : Fit.id,
                            'd' : Expt.id},
                 basis   = ['fit','expt'],
                 constr  = REGEXP(Species['nickname'],BINARY(Fit['dataconst'])))

    pop_fitdata =                                                               \
        Gen(name    = 'pop_fitdata',
            desc    = 'Connect fit to DFT data based on Regex match',
            query   = pfdq,
            tags    = ['fit'],
            actions = [Fit_data(insert   = True,
                                fit      = pfdq['f'],
                                expt     = pfdq['d'])])

    ########################################################################

    def parse_dict(dstr:str,cname:str)->T[str,float]:
        '''Interprets python dictionary with weights, default 1.0'''
        if not dstr:
            return [],[] # type: ignore
        else:
            d = literal_eval(dstr)
            return (cname, d.get(cname,1.0))

    pfc = defaultEnv + Env(Import('ast','literal_eval'))

    pfcq = Query(exprs  = {'f'  : Fit.id,
                           'c'  : Fit['consts'],
                           'cn' : Const['const_name']},
                 basis  = ['Fit','Const'])

    pfcpb = PyBlock(parse_dict,
                    env      = pfc,
                    args     = [pfcq['c'],pfcq['cn']],
                    outnames = ['cn','cw'])

    con = Const(const_name = pfcpb['cn'])

    pop_fitconst =                                                              \
        Gen(name    = 'pop_fitconst',
            desc    = 'Parses weight dictionary',
            query   = pfcq,
            funcs   = [pfcpb],
            tags    = ['fit'],
            actions = [Fit_const(insert       = True,
                                 fit          = pfcq['f'],
                                 const        = con,
                                 const_weight = pfcpb['cw'])])

    pfncq = Query(exprs = {'f' : Fit.id,
                           'c' : Fit['nlconsts'],
                           'n' : Nl_const['nlconst_name']},
                  basis = ['Fit','nonlin_const'])

    pfncpb = PyBlock(parse_dict,
                     env      = pfc,
                     args     = [pfncq[x] for x in 'cn'],
                     outnames = ['n','w'])

    nc = Nl_const(nlconst_name=pfncpb['n'])

    nlinsert = Fit_nl_const(insert          = True,
                            fit             = pfncq['f'],
                            nl_const_weight = pfncpb['w'],
                            nonlin_const    = nc)
    pop_fitnlconst =                                                            \
        Gen(name    = 'pop_fitnlconst',
            desc    = 'Parses weight dictionary' ,
            query   = pfncq,
            funcs   = [pfncpb],
            tags    = ['fit'],
            actions = [nlinsert])

    ########################################################################

    pfs_env = defaultEnv + Env(Import('re','findall'))
    pfsq    = Query(exprs = {'f':Fit.id,'log':Fit['log']})
    pfsns   = ['niter','cost','c_viol']
    pfspb   = PyBlock(parse_fitstep,
                      env      = pfs_env,
                      args     = [pfsq['log']],
                      outnames = pfsns)
    pop_fitstep =                                                      \
        Gen(name    = 'pop_fitstep',
            desc    = 'Analyzes fitting result log, populates Fitstep',
            query   = pfsq,
            funcs   = [pfspb],
            tags    = ['fit', 'parallel'],
            actions = [Fit_step(insert = True,
                                fit    = pfsq['f'],
                                **{x:pfspb[x] for x in pfsns})])

    ########################################################################

    data_dict = dict(
        atomic_contribs      = C('"',Expt['atomic_contribs'],'"'),
        atomic_energies      = C('"',Expt['atomic_energies'],'"'),
        bulk_ratio           = Expt['bulk_ratio'],
        bulk_energy          = Expt['bulk_energy'],
        bulk_contribs        = C('"',Expt['bulk_contribs'],'"'),
        coefs                = C('"',Expt['coefs'],'"'),
        composition          = C('"',Expt['composition'],'"'),
        expt_cohesive_energy = Expt['expt_cohesive_energy'],
        expt_volume          = Expt['expt_volume'],
        expt_bm              = Expt['expt_bm'],
        energy_vector        = Expt['energy_vector'],
        volume_vector        = Expt['volume_vector'],
        contrib_vector       = Expt['contrib_vector'])

    rdq = Query(exprs   = {'f': Fit.id, 'rd':sqldict(data_dict)},
                links   = [Fit_data.r('Fit')],
                basis   = ['fit'],
                aggcols = [Fit['name']])
    raw_data =                                                              \
        Gen(name    = 'fit_rawdata',
            desc    = 'Serialize info stored in expt mapping table',
            query   = rdq,
            tags    = ['fit'],
            actions = [Fit(fit      = rdq['f'],
                           raw_data = rdq['rd'])])

    ########################################################################

    constdict = dict(name   = C('"',Const['const_name'],'"'),
                     s      = COALESCE(Const['s'],'None'),
                     alpha  = COALESCE(Const['alpha'],'None'),
                     val    = COALESCE(Const['val'],'None'),
                     kind   = C('"',Const['kind'],'"'),
                     vec    = C('"',Const['vec'],'"'),
                     weight = Fit_const['const_weight'])


    rcq = Query(exprs    = {'f':Fit.id, 'rc':sqldict(constdict)},
                basis    = ['fit'],
                links    = [Fit_const.r('Fit')],
                opt_attr = [Const['s'],Const['alpha']],
                aggcols  = [Fit['name']])

    raw_const =                                                                 \
        Gen(name    = 'raw_const',
            desc    = 'Serialize info stored in constraint mapping table',
            query   = rcq,
            tags    = ['fit'],
            actions = [Fit(fit       = rcq['f'],
                           raw_const = rcq['rc'])])

    ########################################################################

    nlcdict = dict(name   = C('"',Nl_const['nlconst_name'],'"'),
                   f      = C('"',Nl_const['f'],'"'),
                   df     = C('"',Nl_const['df'],'"'),
                   hess   = C('"',Nl_const['hess'],'"'),
                   ub     = Nl_const['ub'],
                   lb     = Nl_const['lb'],
                   weight = Fit_nl_const['nl_const_weight'])

    rncq = Query(exprs   = {'f' : Fit.id, 'rnc' : sqldict(nlcdict)},
                basis    = ['fit'],
                links    = [Fit_nl_const.r('Fit')],
                aggcols  = [Fit['name']])
    raw_nlconst =                                                           \
        Gen(name    = 'raw_nlconst',
            desc    = 'Serialize info stored in nl_const mapping table',
            query   = rncq,
            tags    = ['fit'],
            actions = [Fit(fit         = rncq['f'],
                           raw_nlconst = rncq['rnc'])])

    ########################################################################

    fit_out = ['log','timestamp','runtime','result','err']

    fit_in  = ['raw_data','raw_const','raw_nlconst','maxiter','basis','bound',
                'initfit','constden','bm_weight','lat_weight']

    fit_env = Env(Import('functionals.fit','Fit')) + defaultEnv

    mfq = Query(exprs     = {'f':Fit.id, **{x:Fit[x] for x in fit_in}},
                opt_attr  = [Fit['raw_data'],
                             Fit['raw_const'],
                             Fit['raw_nlconst']])

    mfpb = PyBlock(fitting,
                   env      = fit_env,
                   args     = [mfq[x] for x in fit_in],
                   outnames = fit_out)
    make_fit =                                                              \
        Gen(name    = 'make_fit',
            desc    = "Scipy minimize: fill out 'Fit['err' iff it fails",
            query   = mfq,
            funcs   = [mfpb],
            tags    = ['fit','parallel'],
            actions = [Fit(fit = mfq['f'],
                           **{x:mfpb[x] for x in fit_out})])

    ########################################################################

    countdict = {'n_const'   : Fit_const,
                 'n_nlconst' : Fit_nl_const,
                 'n_data'    : Fit_data}

    threegens = []

    for x,y in countdict.items():
        q = Query(exprs   = {'f':Fit.id, 'n':COUNT(Literal(1))},
                  aggcols = [Fit.id],
                  basis   = [y.name])

        threegens.append(
            Gen(name    = 'pop_'+x,
                desc    = 'Aggregate and count Fit-related mapping tables',
                query   = q,
                tags    = ['fit'],
                actions = [Fit(fit = q['f'],
                               **{x:q['n']})]))

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
            query   = gq,
            tags    = ['fit'],
            actions = [Globals(insert = True, **{g:gq[g] for g in gcols})])

    ########################################################################
    resid_cols = ['r2_ce','r2_bm','r2_lat','c_viol']

    fm_env = defaultEnv + Env(Import('functionals','FromMatrix'))

    def f(dat:str,con:str,nlcon:str,x:str,n:int)->T[float,float,float,float]:
        return FromMatrix(np.array(loads(x)).reshape((n,n))).costs(dat,con,nlcon)

    rq = Query(exprs = {'f'   : Fit.id,
                        'res' : Fit['result'],
                        'b'   : Fit['basis'],
                        **{g:Globals[g] for g in gcols}},
                basis = ['Fit','Globals'])

    rpb = PyBlock(f,
                  env      = fm_env,
                  args     = [rq[x] for x in gcols+['res','b']],
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

    def bdist(x:str,n:int)->float:
        fx = FromMatrix(np.array(loads(x)).reshape((n,n)))
        r  = np.logspace(-2,2,5).tolist()
        return Functional.diff(fx,BEEF,r,r)

    beef_env = Env(Import('functionals','BEEF','Functional'))

    bdpb = PyBlock(bdist,
                   env  = fm_env + beef_env,
                   args = [fbq['r'],fbq['b']])
    beefdist =                                                                  \
        Gen(name    = 'beefdist',
            desc    = 'Computes Fit.beefdist',
            query   = fbq,
            funcs   = [bdpb],
            tags    = ['fit'],
            actions = [Fit(fit      = fbq['f'],
                           beefdist = bdpb['out'])])

    ########################################################################

    def f_score(r2ce:float, r2bm:float, r2lat:float, cviol:float) -> float:
        return float(2*r2ce + r2bm + r2lat) - log10(cviol)/5

    menv = Env(Import('math','log10'))

    sq   = Query({'f':Fit.id,
                 **{x:Fit[x] for x in resid_cols}})
    spb  = PyBlock(f_score,
                   env  = menv,
                   args = [sq[x] for x in resid_cols])
    score =                                                                     \
        Gen(name    = 'score',
            desc    = 'Computes Fit.score',
            query   = sq,
            funcs   = [spb],
            tags    = ['fit'],
            actions = [Fit(fit   = sq['f'],
                           score = spb['out'])])

    ########################################################################

    def get_lda_viol(x : str, n : int) -> float:
        M = np.array(loads(x)).reshape((n,n))
        return abs(1 - FromMatrix(M).apply(s=0,a=1))

    lvpb = PyBlock(get_lda_viol,
                   env  = fm_env,
                   args = [fbq['r'], fbq['b']])
    lda_viol =                                                                 \
        Gen(name    = 'lda_viol',
            desc    = 'Computes Fit.lda_viol',
            query   = fbq,
            funcs   = [lvpb],
            tags    = ['fit'],
            actions = [Fit(fit=fbq['f'],lda_viol=lvpb['out'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [pop_constraint,pop_nlconstraint,pop_all_fit,pop_fit,pop_fitdata,
            pop_fitconst,pop_fitnlconst,pop_fitstep,raw_data,raw_const,
            raw_nlconst,nconst,nnlconst,ndata,globals,resid,beefdist,score,
            lda_viol, make_fit]

    mod.add(gens)
