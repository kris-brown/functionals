from typing         import Any
from json           import load,dump
from numpy          import array, zeros    # type: ignore
from numpy.linalg   import lstsq           # type: ignore
import sys
from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint # type: ignore

#########################################################################

def main()->None:
    '''Requires three files: params, data, and constraint (.json)'''

    # Read data files
    #-----------------
    with open('params.json','r') as fi:
        params = load(fi)
        pnames = ['n','bound','maxit','normc','ifit']
        n,bound,maxit,normC,ifit = map(params.get,pnames)

    with open('data.json','r') as fi:
        X,target = map(array,load(fi))

    with open('constraint.json','r') as fi:
        c_A,c_lb,c_ub = map(array,load(fi))

    # Initialize
    #-----------
    n2   = n ** 2
    n2m1 = n2 - 1
    zero = zeros(n2)

    options   = {'verbose'                  : 3,
                 'disp'                     : True,
                 'maxiter'                  : maxit,
                 'xtol'                     : 1e-10,
                 'barrier_tol'              : 1e-10,
                 'gtol'                     : 1e-10,
                 'initial_tr_radius'        : 1.0,
                 'initial_constr_penalty'   : 1.0,
                 'initial_barrier_parameter': 0.1}

    # Define inputs for minimize()
    #------------------------------
    bounds    = Bounds([0.] + [-bound] * n2m1, [2.] + [bound] * n2m1)
    cons      = LinearConstraint(c_A,c_lb,c_ub,keep_feasible=False)

    # Make nonlinear constraint (THIS WILL BECOME 'PARAMETERIZED' LATER (HOW?)
    # FOR NOW IT IS ALWAYS (AND ONLY) A REGULARIZATION PENALTY
    #--------------------------
    def norm(x : array) -> float:
        return x @ x
    def norm_grad(x : array) -> array:
        """ d(XiXi)/d(Xj) = 2 Xi d(Xi)/d(Xj) = 2 Xj"""
        return 2 * x

    norm_hess = lambda _,__: zero # hessian is a constant

    norm = NonlinearConstraint(norm, [0], [normC], jac=norm_grad, hess=norm_hess)

    # Generate initial guess
    #------------------------
    x0 = lstsq(X,target,rcond=None)[0] if ifit else zero

    # Objective function and its derivatives
    #---------------------------------------
    def f(x : array) -> float:
        """Objective function for optimization: sum of squares of residual"""
        resid = X @ x - target
        return resid @ resid

    def df(x:array)->array:
        """Jacobian of objective function (derived below in Einstein notation)"""
        residual = X @ x - target
        return 2 * X.T @ residual

    hess = lambda _: zero # hessian is a constant

    # Redirect stdout for the duration of training
    #-----------------
    original   = sys.stdout
    sys.stdout = open('log', 'w')
    res        = minimize(f, x0, method = 'trust-constr', jac = df,
                          constraints = [cons,norm], hess = hess,
                          bounds = bounds, options = options)
    sys.stdout = original

    # Process results
    #-----------------
    outputs = ['fun','execution_time','x','constr_penalty']
    d = {k:getattr(res,k) for k in outputs}
    d['x'] = d['x'].tolist() # can't JSON numpy arrays

    with open('result.json','w') as fi:
        dump(d,fi)


if __name__ == '__main__':
    main()


"""
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
