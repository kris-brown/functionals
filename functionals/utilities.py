# External modules
from typing import List,Tuple,TypeVar,Callable as C, Optional as O
from numpy import zeros,array                        # type: ignore
from numpy.polynomial.legendre import Legendre # type: ignore
import numpy as np  # type: ignore
from scipy.optimize import linprog # type: ignore
# Type Synonyms
A = TypeVar('A')
Binary = C[[float,float],float]

"""
Miscellaneous
"""
################################################################################
def leg(n:int,x:float)->float:
    """
    Value of nth legendre polynomial evaluated at x
    """
    v = zeros(n+1)
    v[n] = 1 # n = 6 ===> [0,0,0,0,0,1]
    return Legendre(v)(x)

def transformed_s(s:float)->float:
    """
    Properties: s=0 -> -1 , s=inf -> 1
    Eq # 2 in mBEEF paper
    """
    assert s >= 0
    kappa = 0.804
    mu    = 10./81
    q     = kappa / mu # Pade approximant to PBEsol Fx(s)
    return 2 * s**2 / (q+s**2) - 1

def transformed_alpha(alpha:float)->float:
    """
    Properties: a=0 -> 1, a = 1 -> 0, a = inf -> -1
    Eq # 3 in mBEEF paper
    MULTIPLIED BY NEGATIVE ONE!!!
    """
    assert alpha >= 0
    return -(1-alpha**2)**3/(1+alpha**3 + alpha**6)

def LegendreProduct(s:float,alpha:float,M:int,N:int)->float:
    """
    s     - reduced density gradient
    alpha - reduced kinetic energy density
    M,N   - L. polynomials
    Eq # 4 in mBEEF paper
    """
    return leg(M,transformed_s(s)) * leg(N,transformed_alpha(alpha))

def coefFunc_(n:int)->C[[array],Binary]:
    """
    Given a certain number of basis Legendre functions to use, we can construct
    a function that takes an array of n^2 coefficients and returns the binary
    function of a functional.
    """
    def f(x:array)->Binary:
        assert len(x)==n**2
        X = x.reshape(n,n)
        def enhancement(s:float,a:float)->float:
            return sum([LegendreProduct(s,a,j,i) for j in range(n) for i in range(n)])
        return enhancement
    return f

def coefFunc(x:array)->Binary:
    """
    Given a certain number of basis Legendre functions to use, we can construct
    a function that takes an array of n^2 coefficients and returns the binary
    function of a functional.
    """
    n = int(len(x)**0.5)
    assert len(x)==n**2
    X = x.reshape(n,n)
    def enhancement(s:float,a:float)->float:
        return sum([LegendreProduct(s,a,j,i) for j in range(n) for i in range(n)])
    return enhancement

########################
def corner(n:int)->List[int]:
    """
    List of indices in a flattened 8x8 matrix ([00,01,02,...07,10,11,...77])
    that correspond to the upper right corner of size n
    """
    return [x for x in range(8**2) if  x % 8 < n and x < 8*n]

def not_corner(n:int)->List[int]:
    """Opposite of corner"""
    return [i for i in range(8**2) if i not in corner(n)]

###################
def flatten(lol: List[List[A]])->List[A]:
    """Convert list of lists to a single list via concatenation"""
    return [item for sublist in lol for item in sublist]



def lst_sq_pulp(A      : np.array,
                   b      : np.array,
                   bounds : List[Tuple[float,float]],
                   A_eq   : O[np.array] = None,
                   b_eq   : O[np.array] = None,
                   A_ub   : O[np.array] = None,
                   b_ub   : O[np.array] = None,
                   maxiter: int   = 1000,
                   tol    : float = .01,
                   verbose: bool  = True
                   ) -> np.array:
    from pulp import LpVariable,LpProblem,LpMinimize,LpSolverDefault,lpDot,value # type: ignore

    n,m  = A.shape

    s1,s2 = (0,),((0,m))

    if A_eq is None or b_eq is None:
        A_eq,b_eq = (np.empty(s2),np.empty(s1))
    if A_ub is None or b_ub is None:
        A_ub,b_ub = (np.empty(s2),np.empty(s1))

    c    = np.concatenate((np.zeros(m),[1])) # objective functional only minimizes v
    A_   = np.c_[A,-np.ones((n,))]           # add a "- v" to each row in Ax (-v) = b
    A_ub = np.c_[A_ub,np.zeros(len(b_ub))]   # make prexisting systems of eqs ignore v
    A_eq = np.c_[A_eq,np.zeros(len(b_eq))]   # make prexisting systems of eqs ignore v
    A_ub = np.vstack((A_ub,A_,-A_))          # add new pseudo-ub eqs to existing ones
    b_ub = np.concatenate((b_ub,b,-b))            # To flip '<' to '>' required * by -1

    bounds += [(0,1000)]                   # add bounds of dummy variable

    b_ub[b_ub == np.inf]=1000

    vars = [LpVariable('x%d'%i,lb,ub) for i,(lb,ub) in enumerate(bounds)]
    prob = LpProblem('myProblem',LpMinimize)
    prob += vars[-1]
    for eq_row,eq_val in zip(A_eq,b_eq):
        prob+= lpDot(eq_row,vars) == eq_val
    for ub_row,ub_val in zip(A_ub,b_ub):
        prob+= lpDot(ub_row,vars)  <= ub_val
    LpSolverDefault.msg = 1

    prob.solve()
    ans = np.array([value(x) for x in vars[:-1]])
    class Dummy(object): x = ans
    return Dummy


##################
def lst_sq_linprog(A      : np.array,
                   b      : np.array,
                   bounds : List[Tuple[float,float]],
                   A_eq   : O[np.array] = None,
                   b_eq   : O[np.array] = None,
                   A_ub   : O[np.array] = None,
                   b_ub   : O[np.array] = None,
                   maxiter: int   = 1000,
                   tol    : float = .01,
                   verbose: bool  = True
                   ) -> np.array:
    """
    Convert an Ax=b least squares problem into a linear programming problem.
    In doing so, we have the ability to impose equality and inequality
        constraints on the solution

    # Define L1 norm cost
    #--------------------
    # To represent |Ax - b|, we need to set:
    # - v < Ax - b < v
    # Where 'v' is our dummy variable, bounded to be a positive number

    # Scipy only wants inequalities of the form A_ub * [x1,x2,...,v] < b_ub
    # So we need two rows for each row in the original Ax = b

    # 1.) Ai1*x1 + Ai2*x2 + ...  - bi < v
    # 2.) Ai1*x1 + Ai2*x2 + ...  - bi > - v

    # The only way to satisfy both of these is if v is zero and Ax - b = 0


    """
    n,m  = A.shape

    s1,s2 = (0,),((0,m))

    if A_eq is None or b_eq is None:
        A_eq,b_eq = (np.empty(s2),np.empty(s1))
    if A_ub is None or b_ub is None:
        A_ub,b_ub = (np.empty(s2),np.empty(s1))

    c    = np.concatenate((np.zeros(m),[1])) # objective functional only minimizes v
    A_   = np.c_[A,-np.ones((n,))]           # add a "- v" to each row in Ax (-v) = b
    A_ub = np.c_[A_ub,np.zeros(len(b_ub))]   # make prexisting systems of eqs ignore v
    A_eq = np.c_[A_eq,np.zeros(len(b_eq))]   # make prexisting systems of eqs ignore v
    A_ub = np.vstack((A_ub,A_,-A_))          # add new pseudo-ub eqs to existing ones
    b_ub = np.concatenate((b_ub,b,-b))            # To flip '<' to '>' required * by -1

    bounds += [(0,np.inf)]                   # add bounds of dummy variable
    options = {'maxiter': maxiter,
               'tol'    : tol,
               'disp'   : verbose,
               'lstsq'  : True,
               'ip'     : True,
               'sym_pos': False,
               'alpha0' : 0.8,
               'beta'   : 0.01}

    # Try to solve linear programming problem, increase tolerance for constraint
    # violation iteratively until we succeed
    complete = False

    b_ub[b_ub == np.inf]=1000
    while not complete:
        try:
            ans = linprog(c=c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=bounds,
                          options=options,method='interior-point')
            assert not isinstance(ans.x,float)
            complete = True
        except AssertionError:
            options['tol']*=2
            print('tol = %f'%options['tol'])
            pass


    ans.x = ans.x[:-1] # discard the dummy value
    return ans
