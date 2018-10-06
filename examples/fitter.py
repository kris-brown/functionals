import numpy as np # type: ignore
from numpy import exp,pi# type: ignore
from scipy.integrate import quad# type: ignore
from scipy.optimize import root# type: ignore


### changes for funtional forms, mu and kappa values ###
def fx(p:float,mu:float,kappa:float,c:float)->float:
    if p>=0.:
        return 1. + kappa - kappa/(1.+(mu*p + c)/kappa)
    else:
        return 1. + kappa  # return here the limit for reduced density gradient -> infinity (or p->infinity)



#initial guess for c
iniguess = 0.1

### changes for funtional forms, mu and kappa values ###

ctmp = [0.0]

def intkernel(r:float)->float:
    mu,kappa = [10./81.,0.29]

    if r<500.:
        p = exp(4.E0/3.E0*r)/(pi*6.)**(2.E0/3.E0)
    else:
        p = -1.

    fxterm = fx(p,mu,kappa,ctmp[0])
    const  = -3.*3.**(1.E0/3.E0)/pi**(2.E0/3.E0)*2.**(1.E0/3.E0)
    return const*exp(-8.E0/3.E0*r)*r**2*fxterm

def fun(c:float)->float:
    ctmp[0] = c
    return -0.3125 - quad(intkernel, 0., np.inf)[0]

sol = root(fun, iniguess)

print(sol)
print('c-value:', sol.x[0])
print('Hydrogen exchange-energy:', fun(sol.x)-0.3125)
print('Hydrogen analytical Hartree energy: 0.3125')
