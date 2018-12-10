# External
from typing    import Dict as D, List as L, Tuple as T, Any, Optional as O
from ast        import literal_eval
from time       import time
from json       import dumps, loads
from io         import StringIO
from contextlib import redirect_stdout
from traceback  import format_exc
from contextlib import contextmanager

import sys

from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint # type: ignore
from numpy          import array,logspace,inf,zeros         # type: ignore
from numpy.linalg   import lstsq                            # type: ignore
import warnings; warnings.filterwarnings("ignore")
from sklearn.linear_model import LinearRegression           # type: ignore

# Internal
from functionals.fit.utilities import LegendreProduct
from dbgen import flatten
from functionals.fit.data import process_data,weight

'''

Try "keep_feasible" = True???
'''
################################################################################



class Fit(object):
    """
    Fit a <n x n> matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """

    def __init__(self,
                 data           : str,
                 constraints    : str,
                 nlconstraints  : str,
                 maxit          : int,
                 n              : int,
                 bound          : float,
                 ifit           : bool,
                 gridden        : int,
                 bm_weight      : float,
                 lat_weight     : float
                 ) -> None:

        self.data = data; self.constraints = constraints;
        self.nlconstraints = nlconstraints; self.maxit = maxit; self.n = n;
        self.bound = bound; self.ifit = ifit;
        self.gridden = gridden;

        self.options = dict( verbose                  = 3,
                             disp                     = True,
                             maxiter                  = maxit,
                             xtol                     = 1e-10,
                             barrier_tol              = 1e-10,
                             gtol                     = 1e-10,
                             initial_tr_radius        = 1.0,
                             initial_constr_penalty   = 1.0,
                             initial_barrier_parameter= 0.1)

        self.c_A,self.c_lb,self.c_ub =                                          \
            self.process_constraints(self.constraints,n,gridden)
        self.nl   = self.process_nlconstraints(self.nlconstraints)
        self.zero = zeros(self.n**2)
        self.X0,self.X1,self.X2,self.Y0,self.Y1,self.Y2 = \
            process_data(self.eval_sql_dict(self.data),n)
        self.X, self.Y = weight(bm_weight,lat_weight,
                                (self.X0,self.X1,self.X2,self.Y0,self.Y1,self.Y2))

    @classmethod
    def process_constraints(cls, cons : str, n : int = 8, gridden : int = 5
                            ) -> T[array,array,array]:
        '''
        Convert linear constraints into a coef matrix + upper/lower bound arrays
        '''
        constraints = cls.eval_sql_dict(cons)
        keys = ['s','alpha','val','kind','vec','weight']
        grid = logspace(-2,2,gridden).tolist()
        coefs,lo,hi = [],[],[] # type: T[list,list,list]

        for s_,a_,val_,kind,vec,weight in [map(c.get,keys) for c in constraints]:
            w   = float(weight)
            assert w > 0
            val = float(val_) # type: ignore
            if   kind == 'eq': lo_,hi_ = val, val
            elif kind == 'lt': lo_,hi_ = -inf, val
            elif kind == 'gt': lo_,hi_ = val, inf
            else: raise ValueError

            if vec:
                loaded = loads(vec)
                lo.append(lo_ * w); hi.append(hi_ * w)
                inds = [x for x in range(8**2) if  x % 8 < n and x < 8*n]
                coefs.append([loaded[i] * w for i in inds])
            else:
                srange = [float(s_)] if s_ else grid
                arange = [float(a_)] if a_ else grid
                card   = len(srange) * len(arange)
                w_     = w / card # reduce weight if many data points
                for s in srange:
                    for a in arange:
                        lo.append(lo_* w_)
                        hi.append(hi_* w_)
                        coefs.append(flatten([[LegendreProduct(s,a,i,j) * w_
                                                for j in range(n)]
                                                    for i in range(n)]))
        c_A, c_lb, c_ub =  map(array,[coefs, lo, hi])
        return c_A,c_lb,c_ub

    @classmethod
    def process_nlconstraints(cls, nlconstraints : str) -> list:
        nonlinconsts = []
        for d in cls.eval_sql_dict(nlconstraints):
            w = d['weight']
            kwargs = dict(fun  = lambda x: w * eval(d['f']),
                          lb   = w * d['lb'],
                          ub   = w * d['ub'],
                          jac  = lambda x: eval(d['df']),
                          hess = lambda x,y: eval(d['hess']))
            nonlinconsts.append(NonlinearConstraint(**kwargs))
        return nonlinconsts

    def fit_result(self)->T[O[str],float,O[int],O[str],str]:

        #--------------------
        redir = StringIO()
        t     = time()

        # Generate initial guess
        #------------------------
        x0 = lstsq(self.X,self.Y,rcond=None)[0] if self.ifit else self.zero

        try:
            with redirect_stdout(redir):
                res = self.fit(x0)
            dt = int(res.execution_time); x = dumps(res.x.tolist())
            return (redir.getvalue(), t, dt, x, '')

        except:
            return None,t,None,None,format_exc()

    def fit(self, guess : array) -> Any:
        '''
        Fit a functional
        '''

        n2m1   = self.n**2 - 1
        bounds = Bounds([0.] + [-self.bound] * n2m1, [2.] + [self.bound] * n2m1)
        no_cons = [LinearConstraint(zeros((1,self.n**2)), zeros(1),zeros(1))]
        cons   = [LinearConstraint(self.c_A, self.c_lb, self.c_ub, keep_feasible = False)] \
                        if len(self.c_lb) else no_cons


        return minimize(self.f, guess, method = 'trust-constr', jac = self.df,
                        constraints = cons + self.nl, hess = self.hess,
                        bounds = bounds, options = self.options)

    def constr_violation(self, x : array) -> float:

        @contextmanager # type: ignore
        def nostdout()->None:

            class DummyFile(object):
                def write(self, x:str)->None: pass
                def flush(self)->None: pass

            save_stdout = sys.stdout
            sys.stdout = DummyFile() #type: ignore
            yield
            sys.stdout = save_stdout

        with nostdout():
            res = self.fit(x)

        return res.constr_violation

    def resid(self, x : array, kind : str) -> array:
        ''' Return Residual '''
        assert kind in ['ce','bm','lat']
        d = dict(ce=(self.X0,self.Y0),bm=(self.X1,self.Y1),lat=(self.X2,self.Y2))
        X,Y = d[kind]
        return X@x - Y

    def r2(self, x : array, kind : str) -> float:
        '''Return R2'''
        assert kind in ['ce','bm','lat']
        d = dict(ce=(self.X0,self.Y0),bm=(self.X1,self.Y1),lat=(self.X2,self.Y2))
        X,Y =  d[kind]
        a,b = (X@x).reshape((-1,1)),Y.reshape((-1,1))

        reg = LinearRegression().fit(a, b)
        return reg.score(a, b)

    # Objective function and its derivatives
    #---------------------------------------
    def f(self, x : array) -> float:
        """Objective function for optimization: sum of squares of residual"""
        resid = self.X @ x - self.Y
        return resid @ resid

    def df(self, x:array)->array:
        """Jacobian of objective function (derived below in Einstein notation)"""
        residual = self.X @ x - self.Y
        return 2 * self.X.T @ residual

    def hess(self,_:array) -> array: return self.zero # hessian is a constant

    @staticmethod
    def eval_sql_dict(x : str) -> list:
        def invert_dict(x:D[str,L]) -> L[D[str,Any]]:
            '''A dict of lists turned into a list of dicts'''
            tup_to_dict = lambda tup: dict(zip(x.keys(), tup))
            tuples      = zip(*x.values())
            return list(map(tup_to_dict,tuples))
        return invert_dict(literal_eval(x.replace('\n',''))) if x else []

"""
BONUS
-----

Derivation of Jacobian in Einstein notation:

        Capital letters = vectors, matrices in [brackets]

        Fi  = residual                                                      --- (DIM: m x 1)
        Aij = basis function coefficients (cols) for each data point (row)  --- (DIM: m x n)
        Rj  = vector of weights corresponding to basis functions            --- (DIM: n x 1)
        Yi  = vector of targets for dot product of rows in Aij and Rj       --- (DIM: m x 1)
        δij = Kroenecker delta

        let: Fj = [A]ji * Ri - Yj
        d(Fj)/d(Rk) = [A]ji * d(Ri)/d(Rk) = [A]ji * δik = [A]jk

        d(FjFj)/d(Rk) = 2 * d(Fj)/d(Rk) * Fj
                      = 2 * [A]jk * Fj

Derivation of constant hessian:
        To see this, take derivative of jac result w/r/t some Rm:

        d(FjFj)/d(Rk)d(Rm) = 2 * [A]jk * d(Fj)/d(Rm) =  2 * [A]jk * [A]jm
"""
