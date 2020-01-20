# External
from typing import Any, Callable as C, Tuple as T, List as L
import numpy as np
from scipy.integrate import quad
from json import dumps
# Internal
from functionals.fit.math import LegProduct
# Type Synonyms
Unary = C[[float], float]
#############################################################################


def fx(i: int, j: int) -> Unary:
    '''input float in callable is reduced gradient'''
    return lambda s: float(LegProduct(s=s, a=0.)[i, j])


# Shorthand for precise thirds
# ----------------------------
_1 = 1.E0/3.E0
_2 = 2.E0/3.E0
_4 = 4.E0/3.E0
_8 = 8.E0/3.E0


def intkernel(i: int, j: int) -> Unary:
    fx_ij = fx(i, j)

    def f(r: float) -> float:
        """
        Return the exchange energy density at a certain distance from nucleus
        of hydrogen atom --- p = s^2
        """
        p = float(np.exp(_4 * min(r, 500))/(np.pi*6.)**(_2))
        s = p**0.5
        pi2 = float(np.pi**(_2))

        # p is unbounded, so there will be a problem when integrating to
        # infinity unless we do something
        return -3.*3.**(_1)/pi2*2.**(_1)*float(np.exp(-r * _8))*r**2 * fx_ij(s)
    return f


def one_time_function() -> L[float]:
    return [quad(intkernel(i, j), 0., np.inf)[0]
            for j in range(8) for i in range(8)]


def h_norm_const() -> T[str, float]:
    '''Returns a 64-element vector which can be dotted with a BEEF coef vector,
       the constraint being equality with -0.3125'''
    return dumps(one_time_function()), -0.3125


def weighted_s() -> None:

    from plotly.graph_objs import Figure, Histogram as Hist
    from numpy.random import choice
    from numpy import array, logspace
    from plotly.offline import plot

    def rpt(x: Any) -> L[Any]:
        return [x() for _ in range(10000)]

    rs = logspace(-2, 1, 100)
    den = array([r**2 * np.exp(-r) for r in rs])
    pden = den / sum(den)
    p = array([np.exp(_4 * min(r, 500))/(np.pi*6.)**(_2) for r in rs])
    s = p**(0.5)
    dist = rpt(lambda: choice(a=s, p=pden))

    data = [Hist(x=dist,
                 nbinsx=10000,
                 name='hydrogen',
                 histnorm='probability')]

    plot(Figure(data=data))


if __name__ == '__main__':
    weighted_s()
