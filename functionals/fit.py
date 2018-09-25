# External
from typing      import Tuple, Dict, List, Callable, Any
import numpy as np                     # type: ignore
from scipy.optimize import minimize,Bounds,BFGS,LinearConstraint        # type: ignore
# Internal Modules
from functionals.functional import FromMatrix,beefcoeff
from functionals.constraint import Constraint
from functionals.data       import Data
from functionals.utilities  import corner

###############################################################################
class Fit(object):
    """
    Data + constraints, sufficient info to solve fitting problem
    """
    def __init__(self,n:int,data:List[Data],constraints:List[Constraint])->None:
        self.n = n; self.constraints = constraints; self.data = Data.merge(data)
        self._res = None; self.weights= None # will get computed only once and stored here
        self.lin = [c for c in self.constraints if c.is_linear]
        self.nonlin = [c for c in self.constraints if not c.is_linear]


    def solve(self)->np.array:
        """
        Find optimal (constrained) coefficients and store in the _res field
        """
        if self._res == None:

            # General solver inputs
            #----------------------
            x0        = [beefcoeff[i] for i in corner(self.n)] # np.array([-9] * self.n ** 2)
            bounds    = Bounds([-2] * self.n ** 2, [2] * self.n ** 2)
            eq_cons   = [c.make(self.n) for c in self.lin]
            ineq_cons = [c.make(self.n) for c in self.nonlin]
            trivial   = [LinearConstraint([0] * self.n ** 2,[0],[0])]
            cons      = (eq_cons+ineq_cons) or trivial
            options   = {'verbose' : 1,'maxiter':1e6,'xtol':1e-14,
                         'barrier_tol':1e-14}

            # unpack objective function and derivatives from Data
            f, df, hess = self.data.obj(), self.data.jac(), self.data.hes()

            self._res = minimize(f, x0, method = 'trust-constr', jac = df,
                                 constraints = cons, hess = hess,
                                 bounds = bounds, options = options)
            self.weights = self._res.x # type: ignore

        return self.weights

    def resid(self)->np.array:
        """ Compute residual """

        return self.data.A(self.n**2) @ self.solve() - self.data.target

    def spaghetti(self)->np.array:
        #Effective number of degrees of freedom gives a scale of the ensemble
        #temperature T something to play with, controls the Spaghetti spread
        #----------------------------------------------------------------------
        Ndeg = 2.5                          # Equation 8 (for actual formula)
        min_cost = 0.5 * np.sum(self.resid()**2)     # Equation 9
        T   = 2. * min_cost / Ndeg          # Equation 11b
        A   = np.dot(self.data.x.T,self.data.x) / T
        l,U = np.linalg.eigh(A)
        L   = np.diag(1./np.sqrt(l))
        M   = np.dot(U,L)  # ensemble information
        return M

    def functional(self,name:str='Fit')->FromMatrix:
        return FromMatrix(self.solve().reshape((self.n,self.n)),name=name)
