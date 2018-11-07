from typing     import Any, Optional as O, Tuple as T, List as L, Dict
from ast        import literal_eval
from io         import StringIO
from time       import time
from json       import dumps
from contextlib import redirect_stdout
from numpy      import array,logspace,inf,vstack,zeros # type: ignore
from numpy.linalg   import lstsq # type: ignore
from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint # type: ignore
import warnings
warnings.filterwarnings("ignore")

# Internal
from functionals.fit.utilities import LegendreProduct
from dbgen                     import flatten

###############################################################################
def fitting(cohesive_data_: str,
            constraints_  : str,
            maxit        : int,
            n            : int,
            bound        : float,
            normC        : float,
            ifit         : bool,
            gridden      : int
           ) -> T[O[str],float,O[int],O[str],O[float],O[float],str]:

    def invert_dict(x:Dict[str,L])->L[Dict[str,Any]]:
        f        = lambda y: dict(zip(x.keys(),y))
        unpacked = zip(*x.values())
        return list(map(f,unpacked))

    ev = lambda x: invert_dict(literal_eval(x.replace('\n','')))

    cohesive_data,constraints = map(ev,[cohesive_data_,constraints_])

    # Convert constraints into matrices
    #------------------------------------
    keys = ['s','alpha','val','kind']
    grid = logspace(-2,2,gridden).tolist()
    coefs,lo,hi = [],[],[] # type: T[list,list,list]


    for s_,a_,val_,kind in [map(c.get,keys) for c in constraints]:
        val    = float(val_) # type: ignore
        srange = [float(s_)] if s_ else grid
        arange = [float(a_)] if a_ else grid
        if   kind == 'eq': lo_,hi_ = val, val
        elif kind == 'lt': lo_,hi_ = -inf, val
        elif kind == 'gt': lo_,hi_ = val, inf
        else: raise ValueError

        for s in srange:
            for a in arange:
                lo.append(lo_); hi.append(hi_)
                coefs.append(flatten([[LegendreProduct(s,a,i,j)
                                        for j in range(n)]
                                            for i in range(n)]))
    c_A, c_lb, c_ub =  map(array,[coefs, lo, hi])


    # Convert "data" input
    #----------------------
    def _process(d:dict)->tuple:
        '''Convert a row data into a row of coefficients and a target'''
        comp  = literal_eval(d['composition']) # type: Dict[int,int]

        def mkarray(xs:list) -> array:
            return array(xs)[:n,:n].reshape((n**2))

        try:

            ac_items        = literal_eval(d['atomic_contribs']).items()
            atomic_contribs = {i:mkarray(xs) for i,xs in ac_items}
            ae_items        = literal_eval(d['atomic_energies']).items()
            atomic_energies = {i:float(e) for i,e in ae_items}

            e_atom = sum([atomic_energies[elem]*num for elem,num in comp.items()])
            x_atom = sum([atomic_contribs[elem]*num for elem,num in comp.items()])

        except: import pdb;pdb.set_trace()

        # allows us to convert xs -> energy
        coefs  = mkarray(array(literal_eval(d['coefs'])).reshape((8,8)))

        rat    = int(d['bulk_ratio'])
        e_bulk = float(d['bulk_energy']) / rat
        x_bulk = mkarray(literal_eval(d['bulk_contribs'])) / rat

        dx     = x_atom - x_bulk # THIS IS WHAT WE ARE FITTING

        # Use the BEEF coefficients from this particular calculator
        ex_atom = coefs @ x_atom
        ex_bulk = coefs @ x_bulk

        # Actual target is JUST the exchange component of formation energy
        nonx_e_atom = e_atom - ex_atom
        nonx_e_bulk = e_bulk - ex_bulk
        expt_co_eng = float(d['target']) # expt E atom - E bulk
        target = expt_co_eng  - (nonx_e_atom - nonx_e_bulk) # should be x_atom - x_bulk

        return (dx,target)

    processed   = map(_process, cohesive_data)
    xs,targets  = zip(*processed)
    X,TARGET    = vstack(xs), array(targets)

    ########
    # MAIN #
    ########
    n2m1        = n**2 - 1
    bounds    = Bounds([0.] + [-bound] * n2m1, [2.] + [bound] * n2m1)
    cons      = LinearConstraint(c_A,c_lb,c_ub,keep_feasible=False)


    zero = zeros(n**2)

    options   = {'verbose'                  : 3,
                 'disp'                     : True,
                 'maxiter'                  : maxit,
                 'xtol'                     : 1e-10,
                 'barrier_tol'              : 1e-10,
                 'gtol'                     : 1e-10,
                 'initial_tr_radius'        : 1.0,
                 'initial_constr_penalty'   : 1.0,
                 'initial_barrier_parameter': 0.1}

    def norm(x : array) -> float:
        return x @ x
    def norm_grad(x : array) -> array:
        """ d(XiXi)/d(Xj) = 2 Xi d(Xi)/d(Xj) = 2 Xj"""
        return 2 * x

    norm_hess = lambda _,__: zero # hessian is a constant

    norm = NonlinearConstraint(norm, [0], [normC], jac=norm_grad, hess=norm_hess)

    # Generate initial guess
    #------------------------

    x0 = lstsq(X,TARGET,rcond=None)[0] if ifit else zero
    # Objective function and its derivatives
    #---------------------------------------
    def f(x : array) -> float:
        """Objective function for optimization: sum of squares of residual"""
        resid = X @ x - TARGET
        return resid @ resid

    def df(x:array)->array:
        """Jacobian of objective function (derived below in Einstein notation)"""
        residual = X @ x - TARGET
        return 2 * X.T @ residual

    hess = lambda _: zero # hessian is a constant


    redir = StringIO()

    try:
        with redirect_stdout(redir):
            res = minimize(f, x0, method = 'trust-constr', jac = df,
                           constraints = [cons,norm], hess = hess,
                           bounds = bounds, options = options)

        # Process results
        #-----------------
        outputs = ['fun','execution_time','x','constr_violation']
        d = {k:getattr(res,k) for k in outputs}
        d['x'] = d['x'].tolist() # can't JSON numpy arrays

        log = redir.getvalue()

        return (log, time(), int(d['execution_time']),
                dumps(d['x']), d['fun'], d['constr_violation'], '')
    except Exception as e:
        return None,time(),None,None,None,None,str(e)



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
