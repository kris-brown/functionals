from typing import List
from numpy import array         # type: ignore
import matplotlib               # type: ignore
import numpy as np              # type: ignore
import matplotlib.pyplot as plt # type: ignore

################################################################################

# Manage fonts
#-------------
font = {'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)
matplotlib.rc('legend', fontsize=14)

# Auxillary functions
#-------------------
def legendreN(x:float,n:int)->List[float]:
    f = np.empty(n, np.float)
    if n>1:
        f[0] = 1.
        f[1] = x
        for i in range(2,n):
            f[i] = 2.*x*f[i-1] - f[i-2] - (x*f[i-1] - f[i-2])/i
    elif n==1:
        f[0] = 1.
    return f

def new_s(s:array,coeffs:array)->array:
    #alpha =
    q = 6.5124              # Pade approximant to PBEsol Fx(s)
    t = 2.*s**2/(q+s**2)-1. # in mBEEF paper (between EQ 5 and 6) --- takes reduced gradient s to a t variable that goes from -1 to 1
    n = len(coeffs)
    return np.array([np.dot(legendreN(x,n),coeffs) for x in t])

# Main Function
#-------------
def plot(M:array,p0:array) -> None:
    s        = np.arange(0.,10.01,0.01)
    mu       = 0.2195149727645171
    kappa    = 0.804

    #PBE
    fxpbe    = 1. + kappa*(1. - 1./(1. + mu*s**2 / kappa))
    # RPBE
    fxrpbe   = 1. + kappa*(1. - np.exp(-mu*s**2 / kappa))

    # PBEsol
    mu       = 10./81.
    kappa    = 0.804
    fxpbesol = 1. + kappa*(1.- 1./(1.+mu*s**2/kappa))

    # SCAN
    # TODO

    # Compute the F curve using the exchange coefficients from fitting
    fxfit = new_s(s,p0)

    # Generate ensemble F vs s plots
    n = np.random.normal(size=p0.shape[0])
    q = np.dot(M,n)+p0
    plt.plot(s,new_s(s,q),color='orange', label='Ensemble')

    for i in range(74):
        n = np.random.normal(size=p0.shape[0])
        q = np.dot(M,n)+p0
        plt.plot(s,new_s(s,q),color='orange',alpha=0.1)

    plt.plot(s, fxrpbe,   'c--', label='RPBE')
    plt.plot(s, fxpbe,    'b',   label='PBE')
    plt.plot(s, fxfit,    'r',   label='Fit to expt solid cohesive energies')
    plt.plot(s, fxpbesol, 'm',   label='PBEsol')

    plt.xlabel('Reduced density gradient')
    plt.ylabel('Exchange enhancement')
    plt.legend(loc='best')
    plt.ylim((0.9,2.0))
    plt.savefig('spaghetti.pdf',bbox_inches='tight')

    print('curvature PBE:',    (2.*fxpbe[1]   - 2.*fxpbe[0])   /s[1]**2)
    print('curvature PBEsol:', (2.*fxpbesol[1]- 2.*fxpbesol[0])/s[1]**2)
    print('curvature fit:',    (2.*fxfit[1]   - 2.*fxfit[0])   /s[1]**2)
