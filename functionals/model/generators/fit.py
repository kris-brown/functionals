# External Modules
from typing import Any, Type, Tuple, List, TYPE_CHECKING

# Internal Modules
if TYPE_CHECKING:
    from dbgen.support.model     import Model

from dbgen.support.get      import (CONST, DESC, INPUT, FUNC, GET)
from dbgen.support.expr     import (AS)
from dbgen.support.funclike import SimpleFunc

from functionals.scripts.parse_fits     import  parse_fits
from functionals.scripts.parse_fitstep  import  parse_fitstep
from functionals.scripts.parse_fitconst import  parse_fitconst

################################################################################

# Description, s, alpha, val, kind
def_consts = dict(
    lda    = ('LDA limit', 0, 1, 1, 'eq'),
    liebox = ('The Lieb Oxford bound', None, None, 1.804, 'lt'),
    scan11 = ('Tighter L.O. bound (1.174) when ⍺=0 (Equation 11 in SCAN paper)',
              None, 0, 1.174, 'lt'),
    pos    = ('Make Fx never negative', None, None, 0, 'gt')
)


def const()->Tuple[List[str],List[str],List[float],List[float],List[float],List[str]]:
    names = list(def_consts.keys())
    descs,ss,alphs,vals,kinds = map(list,zip(*def_consts.values()))
    return names,descs,ss,alphs,vals,kinds # type: ignore

def fit(mod:Type['Model']) -> None:
    # Extract tables
    tabs = ['fit','const','fit_step','fit_const']

    Fit, Const, Fit_step, Fit_const = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        constcols = ['const_name','description','s','alpha','val','kind']
        pop_constraint =                                                        \
            ((Const,*Const.select(constcols[1:])) ==
                GET /FUNC/ SimpleFunc(const, outputs = constcols)
                    /DESC/ 'populate constraints')
        ########################################################################
        outs = ['name','resid','runtime','result','timestamp','basis','norm',
                'initfit','bound','maxiter','constden']

        pop_fit =                                                               \
            ((Fit, Fit.basis, Fit.norm, Fit.initfit, Fit.runtime, Fit.result,
                Fit.maxiter, Fit.constden) ==
                GET /FUNC/ SimpleFunc(parse_fits,outputs=outs) /DESC/ 'Looks in FITPATH')

        ########################################################################
        outs = ['name','niter','cost','c_viol']
        pop_fitstep =                                                               \
            ((Fit_step, Fit_step.cost, Fit_step.c_viol) ==
                GET /INPUT/ [] /FUNC/ SimpleFunc(parse_fitstep,outputs=outs) /DESC/ 'Looks in FITPATH')
        ########################################################################
        outs = ['name','constname']
        pop_fitconst =                                                               \
            ((Fit_const) ==
                GET /INPUT/ [] /FUNC/ SimpleFunc(parse_fitconst,outputs=outs) /DESC/ 'Looks in FITPATH')
        ########################################################################
        ########################################################################
