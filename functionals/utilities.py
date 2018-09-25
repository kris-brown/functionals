from typing import List,TypeVar,Callable as C
from numpy import zeros,array                        # type: ignore
from numpy.polynomial.legendre import Legendre # type: ignore
A = TypeVar('A')
Binary = C[[float,float],float]
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
    kappa = 0.804
    mu    = 10./81
    q     = kappa / mu # Pade approximant to PBEsol Fx(s)
    return 2*s**2 / (q+s**2) - 1

def transformed_alpha(alpha:float)->float:
    """
    Properties: a=0 -> 1, a = 1 -> 0, a = inf -> -1
    Eq # 3 in mBEEF paper
    MULTIPLIED BY NEGATIVE ONE!!!
    """
    return -(1-alpha**2)**3/(1+alpha**3 + alpha**6)

def LegendreProduct(s:float,alpha:float,M:int,N:int)->float:
    """
    s     - reduced density gradient
    alpha - reduced kinetic energy density
    M,N   - L. polynomials
    Eq # 4 in mBEEF paper
    """
    return leg(M,transformed_s(s)) * leg(N,transformed_alpha(alpha))

def coefFunc(n:int)->C[[array],Binary]:
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
