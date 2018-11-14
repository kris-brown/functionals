from typing import Callable
from numpy import exp,pi            # type: ignore
import numpy as np                  # type: ignore
from scipy.integrate import quad    # type: ignore
from scipy.optimize import root     # type: ignore

kappa = 0.804
mu    = 10./81

### changes for funtional forms, mu and kappa values ###
def fx(coefmatrix:np.array,s:float,alpha=0)->float:
    """
    Exchange enhancement at alpha = 0
    (in this case, in the form of MS2...could be changed to be in the form of mBEEF)
    Hartree and Bohr (returns atomic units)
    """
    assert alpha==0
    if s > 1000:
        s = 1000
    if p>=0.:
        return 1. + kappa - kappa/(1.+(mu*p + c)/kappa)
    else:
        return 1. + kappa  # return here the limit for reduced density gradient -> infinity (or p->infinity)

# mu_and_kappa = [10./81.,0.29]
#initial guess for c
iniguess = 0.1
### changes for funtional forms, mu and kappa values ###
# ctmp = [0.0]
from typing import Callable

# input float in callable is reduced gradient
def fx(i:int,j:int)->Callable[[float],float]:
    raise NotImplementedError


def intkernel(i:int,j:int):
    fx_ij = fx(i,j)
    def f(r:float,fx_ij:Callable)->float:
        """
        Return the exchange energy density at a certain distance from nucleus of hydrogen atom

        p = s^2
        """
        if r<500.: p = exp(4.E0/3.E0*r)/(pi*6.)**(2.E0/3.E0)
        else: p = exp(4.E0/3.E0*500.)/(pi*6.)**(2.E0/3.E0)
        s = p**0.5
        #else:      p = -1.
        i,j = i_and_j
        # p is unbounded, so there will be a problem when integrating to infinity unless we do something
        return -3.*3.**(1.E0/3.E0)/pi**(2.E0/3.E0)*2.**(1.E0/3.E0)*exp(-8.E0/3.E0*r)*r**2 * fx_ij(s)
    return f

from numpy import array # type: ignore

def one_time_function()->array:
    return [quad(intkernel(i,j), 0., np.inf)[0]
                    for j in range(8) for i in range(8)]

def fun(c):
    # Target: -0.3125
    # Nonlinear function: quad(intkernel,)
    ctmp[0] = c
    return -0.3125 - quad(intkernel, 0., np.inf)[0]

# sol = root(fun, iniguess)
#
# print(sol)

print( 'c-value:', sol.x[0])
print( 'Hydrogen exchange-energy:', fun(sol.x)-0.3125)
print( 'Hydrogen analytical Hartree energy: 0.3125')

def nonlinear_cost(x:array)->float:
    '''Only have a cost if alpha = 0'''
    return intkernel() + .3125
def nonlinear_grad(x:array)->array:
    '''Gradient is a constant (values differ for Cij)
    return ?
def nonlinear_hess(x:array,y:array)->float:
    return 0
