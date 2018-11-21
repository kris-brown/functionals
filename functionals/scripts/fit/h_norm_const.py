# External
from typing import Callable as C, Tuple as T
from numpy import exp,pi,array,inf  # type: ignore
from scipy.integrate import quad    # type: ignore
from json import dumps
# Internal
from functionals.fit.utilities import LegendreProduct
# Type Synonyms
Unary = C[[float],float]
################################################################################
def fx(i : int, j : int) -> Unary:
    '''input float in callable is reduced gradient'''
    return lambda s: LegendreProduct(s,0,i,j)

# Shorthand for precise thirds
#-----------------------------
_1  = 1.E0/3.E0
_2  = 2.E0/3.E0
_4  = 4.E0/3.E0
_8  = 8.E0/3.E0

def intkernel(i : int, j : int) -> Unary:
    fx_ij = fx(i,j)

    def f(r :float) -> float:
        """
        Return the exchange energy density at a certain distance from nucleus of hydrogen atom
        p = s^2
        """
        p = exp(_4 * min(r,500))/(pi*6.)**(_2)
        s = p**0.5

        # p is unbounded, so there will be a problem when integrating to infinity unless we do something
        return -3.*3.**(_1)/pi**(_2)*2.**(_1)*exp(-r * _8)*r**2 * fx_ij(s)
    return f

def one_time_function()->array:
    return [quad(intkernel(i,j), 0., inf)[0]
                    for j in range(8) for i in range(8)]

def h_norm_const()->T[str,float]:
    '''Returns a 64-element vector which can be dotted with a BEEF coef vector,
       the constraint being equality with -0.3125'''
    return dumps(one_time_function()),-0.3125