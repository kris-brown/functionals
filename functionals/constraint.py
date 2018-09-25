# External Modules
from typing import List,Tuple,Any,Callable as C, Optional as O
from abc import ABCMeta,abstractmethod
from numpy import array # type: ignore
import numpy as np      # type: ignore
from scipy.optimize import LinearConstraint # type: ignore
from scipy.optimize import NonlinearConstraint # type: ignore

# Internal Modules
from functionals.utilities import flatten
from functionals.utilities  import LegendreProduct

# Type synonyms #
Binary = C[[float,float],float]
################################################################################

################
# Helper funcs #
################

def const(x:Any)->C[[float,float],Any]:
    """Take a value and make a binary constant function"""
    return lambda _,__:x

true = const(True)

################
# Main classes #
################
class Constraint(object,metaclass=ABCMeta):
    """
    kind = {'lt','eq','gt','ltgt'}
    """
    @property
    @abstractmethod
    def is_linear(self)->bool: pass

    @abstractmethod
    def make(self,n:int)->Any:
        raise NotImplementedError


    def __init__(self, name : str, desc : str, kind : str = 'eq') -> None:
        assert kind in ['lt','gt','eq','ltgt']
        self.kind = kind; self.name = name; self.desc = desc

    def bounds(self)->Tuple[float,float]:
        v = getattr(self,'val')
        if self.kind == 'eq':
            return v, v
        elif self.kind == 'lt':
            return -np.inf, v
        elif self.kind == 'gt':
            return v, np.inf
        else:
            raise ValueError

class LinConstraint(Constraint):
    """Force the functional to take a value at a particular s,alpha point"""
    @property
    def is_linear(self)->bool: return True

    def __init__(self,
                 name : str, desc : str,
                 s : float, alpha : float, val : float,
                 kind : str = 'eq')->None:
        super().__init__(name,desc,kind)
        self.s = s;self.a = alpha; self.val = val

    def make(self,n:int)->LinearConstraint:
        coefs = [[LegendreProduct(self.s,self.a,i,j) for j in range(n)] for i in range(n)]
        lo,hi = self.bounds()
        return LinearConstraint([flatten(coefs)],[lo],[hi])

class ExplicitConstraint(Constraint):
    def __init__(self,name:str,desc:str,vec:C[[int],array],lb:float,ub:float)->None:
        super().__init__(name,desc,'ltgt')
        self.vec = vec; self.lb =lb ; self.ub= ub
    @property
    def is_linear(self)->bool: return True
    def make(self,n:int)->LinearConstraint:
        return LinearConstraint([self.vec(n)],[self.lb],[self.ub])

cLDA    = LinConstraint('LDA','LDA limit',1,0,1)
#Fx(s,a) < 1.804
cLiebOx = ExplicitConstraint('LiebOx','''
The Lieb Oxford bound

This can be expressed as ∀s,⍺ ∈ ℝ+: Fx(s,⍺) < 1.8

We replace Fx(s,⍺) with: Σx_ij*Bi(s'(s))*Bj(⍺'(⍺))
where "'" represents the transformed quantities (domain [-1,1], rather than ℝ+)

What we want is an expression f(x_ij) < 1.8, but the quantification over all
possible s,⍺ inputs makes this challenging.

However, we observe that for *any*  Legendre polynomial Bi: MAX(Bi(z)) = 1.
Moreover, this maximum occurs at the same input value (+1.0) for every Legendre
Polynomial, so we can say: ∀i,j ∈ ℝ+, y,z ∈ [-1,1]: MAX(Bi(y)Bj(z)) = 1

To express a LT over all possible inputs s,⍺ is more simply expressed as a LT
for the MAX value over all those inputs. Thus:
    MAX(∀s,⍺ ∈ ℝ+: Σx_ij*Bi(s'(s))*Bj(⍺'(⍺))) = MAX(Σx_ij)

We recover a linear constraint, saying the sum of the matrix elements cannot
exceed 1.8

''',vec = lambda n: np.ones(n**2), lb = -np.inf,ub=1.814)

# Fx(s,0) < 1.174

def scan11vec(n:int)->array:
    """
    Given an n x n matrix we want to construct the flattened array corresponding
    to:
    [1 -1 1 -1 ...
     1 -1 1 -1 ...
     ...           ]
    """
    row = [1,-1]*(n//2) + ([1] if n%2 else [])
    return row * n

cSCAN11 = ExplicitConstraint('Scan11','''
Tighter L.O. bound (1.174) when ⍺=0 (Equation 11 in SCAN paper)

This can be expressed as ∀s ∈ ℝ+: Fx(s,0) < 1.17

We replace Fx(s,⍺) with: Σx_ij*Bi(s'(s))*Bj(⍺'(0))

Luckily, ⍺'(0)= -1 lets us simplify: Bj(⍺'(0)) = (-1)^j

To find the worst case scenarios, we really want to constrain the coefficients
so that the MAX of the function (yielded by those coefs) is LT 1.17

MAX(Σx_ij*Bi(s')(-1)^j) < 1.17. The two worst cases are when Bi(s') is ±1, and
the signs of X coefficients cancel out any negative values.

We could express this as a nonlinear constraint now (knowing that s' is
constrained to one of two values), but we can get away with a linear
constraint.

-1.17 < [1,-1,1,-1,...] • x < 1.17

Where the 1's and -1's should change sign whenever the corresponding ⍺ Legendre
polynomial index changes (special considerations needed for n x n matrices where
n is odd)
''',
        vec = scan11vec,lb=-1.174,ub=1.174)

#######################
# NOT YET IMPLEMENTED #
#######################
#
# class NonLinConstraint(Constraint):
#     """
#     Set Fx(s,a) equal/lt/gt to some constant (but not only at a particular point)
#     If the relation holds over some subset of the full domain, specify with a
#     (ℝ+,ℝ+)->Bool function.
#     """
#     @property
#     def is_linear(self)->bool: return False
#
#     def __init__(self,
#                  name   : str, desc : str, val : float,
#                  kind   : str       = 'eq',
#                  dom    : O[Binary] = None
#                 ) ->None:
#         super().__init__(name,desc,kind)
#         self.dom = dom or true; self.val = val
#
#     def make(self,n:int)->NonlinearConstraint:
#
#         def f(x:array)->float:
#             assert len(x) == n**2
#             raise NotImplementedError
#
#         def jac(x:array)->array:
#             """Method of computing the Jacobian matrix (an 1-by-n matrix, where
#             element (1, j) is the partial derivative of f with respect to
#             x[j]).  A callable must have the following signature:
#             jac(x) -> {ndarray, sparse matrix}, shape (m, n). """
#
#             raise NotImplementedError
#
#         def hess(x:array)->array:
#             raise NotImplementedError
#
#         lo,hi = self.bounds()
#
#         return NonlinearConstraint(f, lo, hi, jac=jac, hess=hess)

# (at large s, say s=20): Fx(2s,0)/Fx(s,0) = 0.707
# but, at large s: s'(2s)==s'(s)
# Need to express as derivative condition???
#cSCAN13 = NonLinConstraint('Scan13','Scaling at large s (⍺=0) proportional to 1/sqrt(s) ')

# cGD     = [Constraint('GD',1e10,0,1.174)]#,LinConstraint(1e10,1e10,1.174)]
# cMS2    = Constraint('MS2',1e10,2,1.5)
