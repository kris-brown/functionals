# External Modules
from typing import Any, Type, Tuple as T, List as L
from csv import DictReader
from os import environ
# Internal Modules
from dbgen import (Model, CONST, DESC, INPUT, FUNC, GET, IO, AGG, CONSTS, TAGS,
                    LINKS, OPTION, AFTER, Expr,
                    GROUP_CONCAT as GC, CONCAT as C, AS, COALESCE,
                   SimpleFunc,flatten)

# from functionals.scripts.parse_fits     import  parse_fits
from functionals.scripts.parse_fitstep  import  parse_fitstep
from functionals.scripts.query_fitconst import  query_fitconst
from functionals.scripts.query_fitdata  import  query_fitdata
from functionals.scripts.fitting        import  fitting

################################################################################
################################################################################
################################################################################

functionals_db = '/Users/ksb/Documents/JSON/functionals.json'
################################################################################

# Description, s, alpha, val, kind
def_consts = dict(
    lda    = ('LDA limit', 0, 1, 1, 'eq'),
    liebox = ('The Lieb Oxford bound', None, None, 1.804, 'lt'),
    scan11 = ('Tighter L.O. bound (1.174) when ⍺=0 (Equation 11 in SCAN paper)',
              None, 0, 1.174, 'lt'),
    pos    = ('Make Fx never negative', None, None, 0, 'gt')
)

def const()->T[L[str],L[str],L[float],L[float],L[float],L[str]]:
    names = list(def_consts.keys())
    descs,ss,alphs,vals,kinds = map(list,zip(*def_consts.values()))
    return names,descs,ss,alphs,vals,kinds # type: ignore

################################################################################
fitcols = ['name','constconst','dataconst','basis','norm','initfit','bound',
           'maxiter','constden']
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
def fit(mod:Type['Model']) -> None:
    # Extract tables
    tabs = ['fit','const','fit_step','fit_const','fit_data','cohesive_data']

    Fit, Const, Fit_step, Fit_const, Fit_data, CD = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        constcols = ['const_name','description','s','alpha','val','kind']
        pop_constraint =                                                        \
            ((Const,*Const.select(constcols[1:])) ==
                GET /FUNC/ SimpleFunc(const, outputs = constcols)
                    /DESC/ 'populate constraints')
        ########################################################################
        pop_fit =                                                               \
            ((Fit, *Fit.select(fitcols[1:])) ==
                GET /FUNC/ SimpleFunc(parse_fits,outputs=fitcols)
                    /DESC/ 'Looks in FITPATH for fitting specifications'
                    /IO/   True
                    /TAGS/ 'fit')

        ########################################################################
        pop_fitconst =                                                               \
            (Fit_const ==
                GET /INPUT/  Fit.constconst
                    /FUNC/   SimpleFunc(query_fitconst,
                                        inputs=['pth','constconst'],
                                        outputs=['const_name'])
                    /CONSTS/ {'pth':functionals_db}
                    /DESC/   ''
                    /TAGS/ 'fit')
        ########################################################################
        outs = ['composition','symmetry', 'beef', 'coefs', 'xc', 'pw', 'econv','n_atoms']
        pop_fitdata =                                                               \
            (Fit_data ==
                GET /INPUT/  Fit.dataconst
                    /FUNC/   SimpleFunc(query_fitdata,
                                        inputs=['pth','dataconst'],
                                        outputs=outs)
                    /CONSTS/ {'pth':functionals_db}
                    /DESC/   ''
                    /AFTER/  'cohesive_data' # hack
                    /TAGS/   'fit')

        ########################################################################
        pop_fitstep =                                                               \
            ((Fit_step, Fit_step.cost, Fit_step.c_viol) ==
                GET /INPUT/ Fit.log
                    /FUNC/ SimpleFunc(parse_fitstep,
                                      outputs=['niter','cost','c_viol'])
                    /DESC/ 'Analyzes fitting result log, populates Fitstep'
                    /TAGS/ 'fit')
        ########################################################################


        codata_dict = dict(atomic_contribs = C('"',CD.atomic_contribs,'"'),
                           atomic_energies = C('"',CD.atomic_energies,'"'),
                           bulk_ratio      = CD.bulk_ratio,
                           bulk_energy     = CD.bulk_energy,
                           bulk_contribs   = C('"',CD.bulk_contribs,'"'),
                           coefs           = C('"',CD.coefs,'"'),
                           composition     = C('"',CD.composition,'"'),
                           target          = CD.target)

        raw_data =                                                              \
            (Fit.raw_data == GET /INPUT/ (sqldict(codata_dict) |AS| 'raw_data')
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

        fit_out = ['log','timestamp','runtime','result','resid','c_viol','err']
        fit_in  = ['raw_data','raw_const','maxiter','basis','bound',
                    'norm','initfit','constden']

        make_fit =                                                              \
            (Fit.select(fit_out) ==
                GET /INPUT/ Fit.select(fit_in)
                    /FUNC/ SimpleFunc(fitting,
                                      inputs  = fit_in,
                                      outputs = fit_out)
                    /TAGS/ ['fit','parallel'])
        ########################################################################

        ########################################################################


#
# codata    = C('{',' "atomic_contribs" : [',
#                         GC(C('"',CD.atomic_contribs,'"')),
#                     '], "atomic_energies" : [',
#                         GC(C('"',CD.atomic_energies,'"')),
#                     '], "bulk_ratio" : [',
#                         GC(CD.bulk_ratio),
#                     '], "bulk_energy" : [',
#                         GC(CD.bulk_energy),
#                     '], "bulk_contribs" : [',
#                         GC(C('"',CD.bulk_contribs,'"')),
#                     '], "coefs" : [',
#                         GC(C('"',CD.coefs,'"')),
#                     '], "composition" : [',
#                         GC(C('"',CD.composition,'"')),
#                     '], "target" : [',
#                         GC(CD.target),
#                     ']}') |AS| 'raw_data'
        # constdata = C('{',' "name" : [',
        #                     GC(C('"',Const.const_name,'"')),
        #                 '], "s" : [',
        #                     GC(COALESCE(Const.s,'None')),
        #                 '], "alpha" : [',
        #                     GC(COALESCE(Const.alpha,'None')),
        #                 '], "val" : [',
        #                     GC(COALESCE(Const.val,'None')),
        #                 '], "kind" : [',
        #                     GC(C('"',Const.kind,'"')),']}') |AS| 'raw_const'
