# External Modules
from typing import Type, Tuple as T, List as L
from csv    import DictReader
from math   import log10
from os     import environ
from json   import loads
from numpy  import array,logspace # type: ignore

# Internal Modules
from dbgen import (Model, CONST, DESC, INPUT, FUNC, GET, IO, AGG, CONSTS, TAGS,
                    LINKS, OPTION, BASIS, Expr, Literal, REGEXP, BINARY,
                    GROUP_CONCAT as GC, CONCAT as C, AS, COALESCE, COUNT,
                    SimpleFunc,flatten)

from functionals.fit import FromMatrix,Functional,BEEF

from functionals.scripts.load.parse_fitstep   import  parse_fitstep
from functionals.scripts.fit.fitting          import  fitting
from functionals.scripts.fit.h_norm_const     import  h_norm_const

################################################################################
################################################################################
################################################################################

functionals_db = '/Users/ksb/Documents/JSON/functionals.json'
################################################################################
h_norm_vec,h_norm_val = h_norm_const()

# Description, s, alpha, val, kind, vec
def_consts = dict(
    lda    = ('LDA limit', 0, 1, 1, 'eq', None),
    liebox = ('The Lieb Oxford bound', None, None, 1.804, 'lt',None),
    scan11 = ('Tighter L.O. bound (1.174) when ⍺=0 (Equation 11 in SCAN paper)',
              None, 0, 1.174, 'lt',None),
    pos    = ('Make Fx never negative', None, None, 0, 'gt',None),
    hnorm  = ('Integral satisfies known H solution',
              None,None,h_norm_val,'eq',h_norm_vec)
)

def_nlconsts = dict(
 norm = ('Regularize length of vector', 'x @ x', '2*x', 'zeros(len(x))', 0, 0)
)
################################################################################

def const()->T[L[str],L[str],L[float],L[float],L[float],L[str],L[str]]:
    names = list(def_consts.keys())
    descs,ss,alphs,vals,kinds,vecs = map(list,zip(*def_consts.values()))
    return names,descs,ss,alphs,vals,kinds,vecs # type: ignore

def nlconst()->T[L[str],L[str],L[str],L[str],L[float],L[float]]:
    names = list(def_nlconsts.keys())
    descs,fs,dfs,hess,lbs,ubs = map(list,zip(*def_nlconsts.values()))
    return names,descs,fs,dfs,hess,lbs,ubs # type: ignore

################################################################################
fitcols = ['name','constconst','nlconstconst','dataconst','basis','initfit','bound',
           'maxiter','constden','bm_weight','lat_weight']

def parse_fits()->T[L[str]]:
    vals = {c:[] for c in fitcols} # type: dict
    with open(environ['FITPATH'],'r') as f:
        reader = DictReader(f)
        for row in reader:
            for k,v in vals.items(): v.append(row[k])
    return tuple(map(list,vals.values())) # type: ignore

################################################################################
def sqldict(d:dict)->Expr:
    '''Turn a Python dict into a SQL (serialized) dict'''
    args = [[',',' "%s" : ['%k,GC(v),']'] for k,v in d.items()]
    return C('{',*flatten(args)[1:],'}')

################################################################################
################################################################################
################################################################################

def fit(mod:Type[Model]) -> None:
    # Extract tables
    tabs = ['fit','const','fit_step','fit_const','fit_data','dft_data',
            'nonlin_const','fit_nonlin_const','globals','expt','species']

    Fit, Const, Fit_step, Fit_const, Fit_data, D, Nl_const,Fit_nl_const,\
    Globals,Expt,Species = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        constcols = ['const_name','description','s','alpha','val','kind','vec']
        pop_constraint =                                                        \
            ((Const,*Const.select(constcols[1:])) ==
                GET /FUNC/ SimpleFunc(const, outputs = constcols)
                    /DESC/ 'populate constraints'
                    /TAGS/  'fit')

        ########################################################################
        nlconstcols = ['nlconst_name','description','f','df','hess','lb','ub']
        pop_nlconstraint =                                                        \
            ((Nl_const,*Nl_const.select(nlconstcols[1:])) ==
                GET /FUNC/ SimpleFunc(nlconst, outputs = nlconstcols)
                    /DESC/ 'populate nonlinear constraints'
                   /TAGS/  'fit')

        ########################################################################
        all_fit = dict(name='all',constconst='1',nlconstconst='1',dataconst='1',
                       basis=8,initfit=1,bound=1,maxiter=1000,constden=3,
                       lat_weight=1, bm_weight=1)
        pop_all_fit =                                                           \
            ((Fit, *Fit.select(fitcols[1:])) ==
                GET /CONSTS/ all_fit /TAGS/ 'fit')
        ########################################################################
        pop_fit =                                                               \
            ((Fit, *Fit.select(fitcols[1:])) ==
                GET /FUNC/ SimpleFunc(parse_fits,outputs=fitcols)
                    /DESC/ 'Looks in FITPATH for fitting specifications'
                    /IO/   True
                #    /AFTER/ 'dft_data' # hack
                    /TAGS/ 'fit')
        ########################################################################

                 # Pop'd Table # FK to here # Regexp on this col # Regex pattern
        queries = [(Fit_const,    Const,    Const.const_name,      'constconst'),
                   (Fit_nl_const, Nl_const, Nl_const.nlconst_name, 'nlconstconst'),
                   (Fit_data,     D,        Species.nickname,      'dataconst')]

        pop_fitconst,pop_fitnlconst,pop_fitdata =   [
            (tab1 ==
                GET /BASIS/  [Fit,tab2]
                    /CONST/  REGEXP(regcol,BINARY(Fit._get(pat)))
                    /DESC/   'Applies a fit-specified query to the DB'
                    /TAGS/ 'fit')
            for tab1,tab2,regcol,pat in queries]

        ########################################################################
        pop_fitstep =                                                               \
            ((Fit_step, Fit_step.cost, Fit_step.c_viol) ==
                GET /INPUT/ Fit.log
                    /FUNC/ SimpleFunc(parse_fitstep,
                                      outputs=['niter','cost','c_viol'])
                    /DESC/ 'Analyzes fitting result log, populates Fitstep'
                    /TAGS/ ['fit', 'parallel'])
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
            contrib_vector       = D.contrib_vector
                           )

        raw_data =                                                              \
            (Fit.raw_data == GET /INPUT/ (sqldict(data_dict) |AS| 'raw_data')
                                 /LINKS/ Fit_data.Fit
                                 /AGG/   Fit
                                 /TAGS/ 'fit')

        ########################################################################
        constdict = dict(name   = C('"',Const.const_name,'"'),
                         s      = COALESCE(Const.s,'None'),
                         alpha  = COALESCE(Const.alpha,'None'),
                         val    = COALESCE(Const.val,'None'),
                         kind   = C('"',Const.kind,'"'))

        raw_const =                                                             \
            (Fit.raw_const == GET /INPUT/  (sqldict(constdict) |AS| 'raw_const')
                                  /LINKS/  Fit_const.Fit
                                  /AGG/    Fit
                                  /OPTION/ [Const.s,Const.alpha]
                                  /TAGS/   'fit')
        ########################################################################
        nlcdict = dict(name = C('"',Nl_const.nlconst_name,'"'),
                       f    = C('"',Nl_const.f,'"'),
                       df   = C('"',Nl_const.df,'"'),
                       hess = C('"',Nl_const.hess,'"'),
                       ub   = Nl_const.ub,
                       lb   = Nl_const.lb)

        raw_nlconst =                                                           \
            (Fit.raw_nlconst ==
                GET /INPUT/  (sqldict(nlcdict) |AS| 'raw_nlconst')
                    /LINKS/  Fit_nl_const.Fit
                    /AGG/    Fit
                    /TAGS/   'fit')

        ########################################################################

        fit_out = ['log','timestamp','runtime','result','err']

        fit_in  = ['raw_data','raw_const','raw_nlconst','maxiter','basis','bound',
                    'initfit','constden','bm_weight','lat_weight']

        make_fit =                                                              \
            (Fit.select(fit_out) ==
                GET /INPUT/ Fit.select(fit_in)
                    /FUNC/  SimpleFunc(fitting,
                                       inputs  = fit_in,
                                       outputs = fit_out)
                    /OPTION/ [Fit.raw_data, Fit.raw_const,Fit.raw_nlconst]
                    /TAGS/   ['fit','parallel'])

        ########################################################################

        countdict = {Fit.n_const   : Fit_const,
                     Fit.n_nlconst : Fit_nl_const,
                     Fit.n_data    : Fit_data}

        nconst,nnlconst,ndata =  [                                                            \
            (x == GET /INPUT/ COUNT(Literal(1))
                      /BASIS/ y
                      /AGG/   Fit
                      /DESC/  'Count'
                      /TAGS/  'fit')
                for x,y in countdict.items()]
        ########################################################################
        gcols = ['all_data','all_constraints','all_nlconstraints']
        fcols = ['raw_data','raw_const','raw_nlconst']
        zipped = zip(fcols,gcols)
        globals =                                                               \
            ((Globals,Globals.select(gcols)) ==
                GET /INPUT/ [Fit._get(x)|AS|y for x,y in zipped]
                    /CONST/ (Fit.name == 'default')
                    /TAGS/  'fit')

        ########################################################################
        resid_cols = ['r2_ce','r2_bm','r2_lat','c_viol']

        def f(dat:str,con:str,nlcon:str,x:str,n:int)->T[float,float,float,float]:
            return FromMatrix(array(loads(x)).reshape((n,n))).costs(dat,con,nlcon)

        resid =                                                                 \
            (Fit.select(resid_cols) == GET /INPUT/ (Globals.select(gcols)
                                      +[Fit.result,Fit.basis])
                              /BASIS/ [Fit,Globals]
                              /TAGS/ 'fit'
                              /DESC/ 'Use all available data and constraints to evaluate functional'
                              /FUNC/ SimpleFunc(f,gcols+['result','basis'],
                                                 outputs=resid_cols))
        ########################################################################
        def bdist(x:str,n:int)->float:
            fx = FromMatrix(array(loads(x)).reshape((n,n)))
            r  = logspace(-2,2,5).tolist()
            return Functional.diff(fx,BEEF,r,r)

        beefdist =                                                              \
            (Fit.beefdist == GET /INPUT/ [Fit.result,Fit.basis]
                                 /FUNC/  bdist
                                 /TAGS/  'fit')

        ########################################################################

        def f_score(r2ce:float,r2bm:float,r2lat:float,cviol:float)->float:

            return float(2*r2ce + r2bm + r2lat) - log10(cviol)
        score =                                                                 \
            (Fit.score == GET /INPUT/ Fit.select(resid_cols)
                              /FUNC/ SimpleFunc(f_score,resid_cols,['score'])
                              /TAGS/  'fit')
