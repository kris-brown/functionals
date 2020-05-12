# External
import numpy as np
from scipy.integrate import quad

#############################################################################


def hydro() -> np.ndarray:
    talpha0 = 1.

    mu = 10./81.
    q = 0.804/mu

    def leg(n: int, x: float) -> np.ndarray:
        if(n < 1):
            return np.array([])
        elif(n == 1):
            return np.array([x])
        r = np.empty(n, np.float)
        r[0] = 1.
        r[1] = x
        for i in range(2, n):
            r[i] = ((2.*i-1.)*x*r[i-1] - (i-1.)*r[i-2]) / i
        return r

    leg8talpha0 = leg(8, talpha0)

    # Hydrogen self-interaction cancellation (alpha=0)
    hydrogen8 = np.empty(8, np.float)
    for i in range(8):
        ind = np.zeros(i+1, np.float)
        ind[-1] = 1.

        def intkernel(r: float) -> float:
            if r < 500.:
                p = np.exp(4.E0 / 3.E0 * r) / (np.pi * 6.) ** (2.E0 / 3.E0)
                t_s = 2. * p / (q + p) - 1.
                return -3. * 3. ** (1.E0 / 3.E0) / np.pi ** (2.E0 / 3.E0) * 2. ** (1.E0 / 3.E0) * np.exp(-8.E0 / 3.E0 * r) * r ** 2 * np.polynomial.legendre.legval(t_s, ind)
            else:
                return -3. * 3. ** (1.E0 / 3.E0) / np.pi ** (2.E0 / 3.E0) * 2. ** (1.E0 / 3.E0) * np.exp(-8.E0 / 3.E0 * r) * r ** 2

        hydrogen8[i] = quad(intkernel, 0., np.inf,
                            epsabs=1e-12, epsrel=1e-12)[0]

    # weights for Hydrogen exchange
    return np.outer(leg8talpha0, hydrogen8).flatten()


# def weighted_s() -> None:

#     from plotly.graph_objs import Figure, Histogram as Hist
#     from numpy.random import choice
#     from numpy import array, logspace
#     from plotly.offline import plot

#     def rpt(x: Any) -> L[Any]:
#         return [x() for _ in range(10000)]

#     rs = logspace(-2, 1, 100)
#     den = array([r**2 * np.exp(-r) for r in rs])
#     pden = den / sum(den)
#     p = array([np.exp(_4 * min(r, 500))/(np.pi*6.)**(_2) for r in rs])
#     s = p**(0.5)
#     dist = rpt(lambda: choice(a=s, p=pden))

#     data = [Hist(x=dist,
#                  nbinsx=10000,
#                  name='hydrogen',
#                  histnorm='probability')]

#     plot(Figure(data=data))


# if __name__ == '__main__':
#     weighted_s()
