# External Modules
from typing import Any, Callable as C, List as L, Dict as D
import abc
import json
import os
import numpy as np
import functools
import plotly.graph_objs as go
import functionals.fit.math as math
import plotly

binary = C[[float, float], float]
bvec = C[[float, float], np.ndarray]
###############################################################################
###############################################################################
# Functional constants
# --------------------
mu_pbe = 0.2195149727645171
kappa = 0.804
mu = 10. / 81

###############################################################################


def flatten(lol: L[Any]) -> L[Any]:
    return [item for sublist in lol for item in sublist]


# Functional classes
# ------------------


class Functional(object, metaclass=abc.ABCMeta):
    """
    Something that implements a function (s,⍺)->Fx :: (ℝ+,ℝ+)->ℝ+
    """
    @property
    @abc.abstractmethod
    def name(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def mgga(self) -> bool:
        pass

    @abc.abstractmethod
    @functools.lru_cache(maxsize=1000, typed=False)
    def apply(self, s: float, a: float) -> float:
        raise NotImplementedError

    def plot(self, color: str, long: bool = False) -> L[Any]:
        ss = np.arange(0., 100 if long else 5, 1 if long else 0.1)
        alphas = np.array([0, 1]) if self.mgga else [1]
        styles = ['solid', 'dot', 'dash', 'dashdot']
        out = []
        for sty, a in zip(styles, alphas):
            lab = self.name + r' (α=%d)' % a if self.mgga else self.name
            out.append(go.Scatter(
                x=ss, y=[self.apply(s, a) for s in ss],
                mode='lines', name=lab,
                line=dict(dash=sty, color=color, shape='spline')))
        return out


# Different ways of creating a functional


class FromFunc(Functional):
    """ Directly specify the enhancement factor formula """

    def __init__(self, name: str, f: binary, mgga: bool = True) -> None:
        self._name = name
        self.f = f
        self._mgga = mgga

    @property
    def name(self) -> str:
        return self._name

    @property
    def mgga(self) -> bool:
        return self._mgga

    @functools.lru_cache(maxsize=100, typed=False)
    def apply(self, s: float, a: float) -> float:
        return self.f(s, a)


class FromMatrix(Functional):
    """
    Implicitly define Fx formula in terms of a 2D matrix corresponding to
    coefficients for Legendre Polynomials
    """

    def __init__(self, A: np.ndarray,  # a1: float = 4.9479, msb: float = 1.,
                 name: str = None) -> None:
        self.A = np.array(A).reshape((8, 8))
        # self.a1 = a1 self.msb = msb
        assert self.A.shape == (8, 8)
        self._name = name or '<no name>'

    @staticmethod
    def frompath(n: str) -> 'FromMatrix':
        fi = 'data/beefs/%s.json' % n
        beefpth = '/' + os.path.join(*__file__.split('/')[:-3], fi)
        with open(beefpth, 'r') as f:
            beefcoeff = np.array(json.load(f))
        return FromMatrix(beefcoeff.reshape(8, 8), name=n)

    @property
    def x(self) -> np.ndarray:
        return self.A.flatten()

    @property
    def name(self) -> str:
        return self._name

    @property
    def mgga(self) -> bool:
        return True

    @functools.lru_cache(maxsize=100, typed=False)
    def apply(self, s: float, a: float,) -> float:
        return float(np.sum(math.LegProduct(s=s, a=a, A=self.A)))

    def plot2d(self, alpha: float = 1., plot: bool = False) -> go.Figure:
        import functionals.fit.constraint as constr
        a = alpha

        def d(s: float) -> float:
            return float(constr.deriv_constr(
                math.t_s(s), math.t_alpha(a)) @ self.x)
        xs = np.linspace(0, 3, 100)
        fs, zs = map(np.array, zip(*[(self.apply(s, a), d(s)) for s in xs]))
        # zs = zs/np.max(np.abs(zs))
        data = [go.Scatter(x=xs, y=fs, name='fx @ α=%f' % alpha),
                go.Scatter(x=xs, y=zs, name='dfx/dŝ @ α=1')]
        fig = go.Figure(data=data)
        if plot:
            plotly.offline.plot(fig)
        return fig

    def plot3d(self, smax: float = 3., amax: float = 3, deriv: str = 'fx',
               plot: bool = False) -> go.Figure:
        import functionals.fit.constraint as constr
        N = 50

        dvs = dict(fx=math.LegProduct, d1=constr.deriv_constr,
                   d2=constr.deriv2_constr, d3=constr.deriv3_constr,
                   ad2=constr.deriv2_alpha)  # type: D[str, bvec]

        def d(s: float, a: float) -> float:
            return float(dvs[deriv](s, a).flatten() @ self.x)

        ss, als, fs, i = zip(*[(s, a, self.apply(s, a), d(s, a))
                               for s in np.linspace(0, smax, N)
                               for a in np.linspace(0, amax, N)])
        lowi, hii = min(i), max(i)
        if lowi > 0:
            colorscale = 'bluered'
        else:
            denom = hii - lowi
            zi = -lowi / denom
            if denom:
                colorscale = [[0., 'blue'], [zi / 2, 'green'],  # type: ignore
                              [zi, 'yellow'], [zi + hii / 2 / denom, 'red'],
                              [1, 'black']]
            else:
                colorscale = [[0., 'white'], [1, 'black']]  # type: ignore
        data = [go.Mesh3d(x=ss, y=als, z=fs, colorscale=colorscale,
                          intensity=i)]
        scene = dict(xaxis=dict(title='s'), yaxis=dict(title='alpha'),
                     zaxis=dict(title='F(s,a)'),
                     )
        title = '%s: %s' % (deriv, self.name)
        fig = go.Figure(data=data, layout=dict(title=title, scene=scene))
        if plot:
            plotly.offline.plot(fig)
        return fig


def plots(ps: L['Functional'], plot: bool = False, long: bool = False
          ) -> go.Figure:
    assert len(ps) < 7
    cs = ['rgb(255,0,0)', 'rgb(0,255,0)', 'rgb(0,0,255)', 'rgb(153,0,153)',
          'rgb(0,0,0)', 'rgb(255,255,0)']
    data = flatten([p.plot(color=c, long=long) for p, c in zip(ps, cs)])
    layout = go.Layout(title='Functionals', xaxis=dict(title='s'),
                       yaxis=dict(title='Fx'))
    fig = go.Figure(data=data, layout=layout)
    if plot:
        plotly.offline.plot(fig)
    return fig

###############################################################################
###############################################################################
###############################################################################
# Functional instances
# --------------------

# Functionals generated by Explicit enhancement factor functions
# --------------------------------------------------------------


PBE = FromFunc('PBE', mgga=False,
               f=lambda s, _: 1. + kappa * (1. - 1. / (1. + mu_pbe * s**2 / kappa)))
RPBE = FromFunc('RPBE', mgga=False,
                f=lambda s, _: float(1. + kappa * (1. - np.exp(-mu_pbe * s**2 / kappa))))
PBEsol = FromFunc('PBEsol', mgga=False,
                  f=lambda s, _: 1. + kappa * (1. - 1. / (1. + mu * s**2 / kappa)))


def fxSCAN(s: float, alpha: float) -> float:
    # Scan-specific constants
    h0x, c1x, c2x, dx, b3, k1, a1 = 1.174, 0.667, 0.8, 1.24, 0.5, 0.065, 4.9479
    b2 = (5913 / 405000)**0.5
    b1 = (511 / 13500) / (2 * b2)
    b4 = mu**2 / k1 - 1606 / 18225 - b1**2
    # Edge conditions with numerical instability
    assert s >= 0 and alpha >= 0
    if s < 0.01:
        s = 0.01
    if abs(alpha - 1) < 0.01:
        alpha = 1.001
    # Intermediate values
    th_1a = float(np.heaviside(1. - alpha, 0.5))
    th_a1 = float(np.heaviside(alpha - 1., 0.5))
    s2 = s**2
    x_ = (b1 * s2 + b2 * (1. - alpha) * np.exp(-b3 * (1 - alpha)**2))**2
    x = mu * s2 * (1. + b4 * s2 / mu) + x_
    h1x = 1. + k1 - k1 / (1. + x / k1)
    gx = 1. - np.exp(-a1 * s**(-0.5))
    fx = np.exp(-c1x * alpha / (1. - alpha)) * th_1a - dx * np.exp(c2x / (1 - alpha)) * th_a1
    # Main output
    return float((h1x + fx * (h0x - h1x)) * gx)


SCAN = FromFunc('SCAN', fxSCAN)


def fxMS2(s: float, alpha: float) -> float:
    k, c, b = 0.504, 0.14601, 4.0
    p = s**2
    F1x = 1 + k - k / (1 + mu * p / k)
    F0x = 1 + k - k / (1 + (mu * p + c) / k)
    f = (1 - alpha**2) / (1 + alpha**3 + b * alpha**6)
    return F1x + f * (F0x - F1x)

# (1 + .504 - .504/(1+(10/81)*x^2/.504)) + ((1-y**2) / (1 + y**3 + 4*y**6))*((1 + .504 - .504/(1+((10/81)*x^2+0.14601)/.504))-(1 + .504 - .504/(1+(10/81)*s**2/.504)))


def ms2_s_curv(s: float, a: float) -> float:
    return (16.7993 - 12.3452 * s**2 / (s**2 + 4.0824)**3
            + ((12.3452 * (s**2 - 1.75503) * (a**2 - 1))
               / ((s**2 + 5.26508)**3 * (4 * a**6 + a**3 + 1))))


MS2 = FromFunc('MS2', fxMS2)

# From data
# ---------
BEEF = FromMatrix.frompath('mbeef')
PBESOL = FromMatrix.frompath('pbesol')


if __name__ == '__main__':
    '''Mbeef at s=0,alpha=1 is 1.03'''
    BEEF.A = np.multiply(
        BEEF.A, [[3**(i + j) / 1500 for j in range(8)] for i in range(8)])
    BEEF.A[0][0] += 1
    plotly.offline.plot(plots([PBE, RPBE, PBEsol, BEEF]))

    # beef = Fsmooth()
    # for x in ['fx']:
    #     fig = beef.plot3d(deriv=x, smax=4)
    #     plotly.offline.plot(fig)
