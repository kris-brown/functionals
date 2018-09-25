# External
from typing import Tuple, Dict, List, Callable, Any, Optional as O
from abc import ABCMeta,abstractmethod
import numpy as np # type: ignore
from scipy.optimize import minimize,Bounds,BFGS,LinearConstraint,NonlinearConstraint,linprog # type: ignore
# Internal Modules
from functionals.functional import FromMatrix,beefcoeff
from functionals.constraint import Constraint,MergedConstraint
from functionals.data       import Data
from functionals.utilities  import corner,flatten,lst_sq_linprog

################################################################################
class Fit(metaclass=ABCMeta):
    """
    Fit a <n x n> matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """

    @abstractmethod
    def solve(self,*args:Any,**kwargs:Any)->np.array:
        raise NotImplementedError

    def __init__(self,n:int,data:List[Data],constraints:List[Constraint],bound:float=0.1)->None:
        self.n = n; self.constraints = constraints; self.bound = bound
        self.data       = Data.merge(data);
        self.constraint = MergedConstraint(constraints)

        # Will get computed only once and stored here
        self._res   = None  # type: O[np.array]
        self.weights = None # type: O[np.array]

    def trivial(self)->List[LinearConstraint]:
        """Make a trivial constraint (Solver fails if there are NO constriants)"""
        return [LinearConstraint([0] * self.n ** 2,[0],[0])]

    def spaghetti(self)->np.array:
        #Effective number of degrees of freedom gives a scale of the ensemble
        #temperature T something to play with, controls the Spaghetti spread
        #----------------------------------------------------------------------
        Ndeg = 2.5                          # Equation 8 (for actual formula)
        min_cost = 0.5 * np.sum(self.functional().resid(self.data)[0]**2)     # Equation 9
        T   = 2. * min_cost / Ndeg          # Equation 11b
        A   = np.dot(self.data.x.T,self.data.x) / T
        l,U = np.linalg.eigh(A)
        L   = np.diag(1./np.sqrt(l))
        M   = np.dot(U,L)  # ensemble information
        return M

    def functional(self,name:str='Fit',maxiter:int=1000,verbose:bool=True)->FromMatrix:
        x = self.solve(maxiter=maxiter,verbose=verbose).reshape((self.n,self.n))
        print('fitted coeffs = ',x)
        return FromMatrix(x,name=name)

class LinFit(Fit):
    """
    Use linear programming (no regularization possible)
    See "Example: 1-norm approximation" in the cvxopt user guide

    """
    def solve(self,maxiter:int=1000,verbose:bool=True)->Tuple[np.array,float]: # type: ignore
        """"""
        if self._res == None:
            # Solve
            #-----
            n2 = self.n**2

            A = self.data.A(n2)
            b = self.data.target
            bounds    = [(0.,2.)]+[(-self.bound,self.bound)]*(n2-1)
            A_eq,b_eq,A_ub,b_ub = self.constraint.linprog_matrices(self.n)

            self._res = lst_sq_linprog(A=A,b=b,bounds=bounds,A_eq=A_eq,b_eq=b_eq,
                                       A_ub=A_ub,b_ub=b_ub,verbose=verbose,
                                       maxiter=maxiter)

            self.weights = self._res.x # type: np.array

            import pdb;pdb.set_trace()
            print('weights',self.weights.reshape(self.n,self.n))
        return self.weights # type: ignore

class NonLinFit(Fit):
    """
    Data + constraints, sufficient info to solve fitting problem
    """
    def __init__(self, n : int, data : List[Data],
                 constraints : List[Constraint],
                 bound: float = .1, norm : float = 0.1,
                 initfit : bool = True) -> None:
        self.norm = norm ; self.initfit = initfit
        super().__init__(n,data,constraints,bound)

    def solve(self, maxiter:int=1000,verbose : bool = True) -> np.array: # type: ignore
        """
        Find optimal (constrained) coefficients and store in the _res field

        initfit - If true, make initial guess the result of the least squares
                  solution
        """
        if self._res == None:

            # General solver inputs
            #----------------------
            if self.initfit:
                x0 = self.data.linreg(self.n)[0] # least square fit, ignore constraints
            else:
                x0 = np.zeros(self.n**2)

            n2m1 = self.n ** 2 - 1

            bounds    = Bounds([0.]+[-self.bound] * n2m1, [2.]+[self.bound] * n2m1)

            if self.constraints:
                cons = MergedConstraint(self.constraints).make(self.n) #flatten([c.make(self.n) for c in self.constraints])
            else:
                cons = self.trivial()

            options   = {'verbose'                  : 3 if verbose else 0,
                         'disp'                     : verbose,
                         'maxiter'                  : maxiter,
                         'xtol'                     : 1e-10,
                         'barrier_tol'              : 1e-10,
                         'gtol'                     : 1e-10,
                         'initial_tr_radius'        : 1.0,
                         'initial_constr_penalty'   : 1.0,
                         'initial_barrier_parameter': 0.1}

            # unpack objective function and derivatives from Data
            f, df, hess = self.data.obj(), self.data.jac(), self.data.hes()
            # functions for normalization
            def norm(x:np.array)->float: return x @ x
            def norm_grad(x:np.array)->np.array:
                """ d(XiXi)/d(Xj) = 2 Xi d(Xi)/d(Xj) = 2 Xj"""
                return 2 * x
            def norm_hess(x:np.array,_:np.array)->np.array:
                return np.zeros(len(x))
            norm = NonlinearConstraint(norm,[0],[self.norm],jac=norm_grad,hess=norm_hess)

            # Calculate
            self._res = minimize(f, x0, method = 'trust-constr', jac = df,
                                 constraints = [cons,norm], hess = hess,
                                 bounds = bounds, options = options)

            self.weights = self._res.x
        return self.weights
