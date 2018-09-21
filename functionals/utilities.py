from numpy import zeros,multiply               # type: ignore
from numpy.polynomial.legendre import Legendre # type: ignore
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
    Properties: s=0 -> -1 , s=inf -> 2
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
