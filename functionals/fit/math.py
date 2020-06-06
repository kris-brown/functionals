import pdb
from typing import Callable as C, List as L
import numpy as np
import scipy.optimize as opt
Unary = C[[float], float]
######################################################
q = 0.804 / (10. / 81.)


# GENERAL DERIVATIVE FUNCTIONS

# def chain_rule(df_dgx: Unary, dg_dx: Unary, gx: float) -> Unary:
#     '''Return function which evaluates d(f(g(x)))/dx at a point x'''
#     def fun(x: float) -> float:
#         return dg_dx(x) * df_dgx(gx)
#     return fun


# def chain_rule2(df2_dgx2: Unary, df_dgx: Unary, dg2_dx2: Unary, dg_dx: Unary,
#                 gx: float) -> Unary:
#     '''evaluate d^2(f(g(x)))/dx^2 at a point
#     https://www.wolframalpha.com/input/?i=second+derivative+of+f%28g%28x%29%29
#     '''
#     def f(x: float) -> float:
#         return dg_dx(x)**2 * df2_dgx2(g(x)) + dg2_dx2(x)*df_dgx(gx)
#     return f


# def chain_rule3(df3_dgx3: Unary, df2_dgx2: Unary, df_dgx: Unary,
#                 dg3_dx3: Unary, dg2_dx2: Unary, dg_dx: Unary, g: Unary,
#                 ) -> Unary:
#     '''evaluate d^3(f(g(x)))/dx^3 at a point
#     https://www.wolframalpha.com/input/?i=third+derivative+of+f%28g%28x%29%29
#     '''
#     def f(x: float) -> float:
#         t1 = df3_dgx3(g(x)) * dg_dx(x)**3
#         t2 = 3*dg_dx(x)*dg2_dx2(x)*df2_dgx2(g(x))
#         t3 = dg3_dx3(x)*df_dgx(g(x))
#         return t1 + t2 + t3
#     return f


# def prod_rule(dfdx: Unary, f: Unary, dgdx: Unary, g: Unary) -> Unary:
#     '''Evaluate d(f(x)g(x))/dx'''
#     def z(x: float) -> float:
#         return dfdx(x)*g(x) + dgdx(x)*f(x)
#     return z


# def prod_rule2(df2dx2: Unary, dfdx: Unary, f: Unary,
#                dg2dx2: Unary, dgdx: Unary, g: Unary) -> Unary:
#     def z(x: float) -> float:
#         return g(x)*df2dx2(x) + 2*dfdx(x)*dgdx(x) + f(x)*dg2dx2(x)
#     return z


# def prod_rule3(df3dx3: Unary, df2dx2: Unary, dfdx: Unary, f: Unary,
#                dg3dx3: Unary, dg2dx2: Unary, dgdx: Unary, g: Unary) -> Unary:
#     def z(x: float) -> float:
#         t1 = df3dx3(x)*g(x)
#         t2 = 3*df2dx2(x)*dgdx(x)
#         t3 = 3*dfdx(x)*dg2dx2(x)
#         t4 = dg3dx3(x)*f(x)
#         return t1+t2+t3+t4
#     return z


# BEEF specific functions

def t_s(s: float) -> float:
    '''Transformation of s [0,inf) to [-1,1]'''
    if np.isinf(s):
        return 1.
    else:
        return 2 * s**2 / (q + s**2) - 1


def t_alpha(alpha: float) -> float:
    '''Transformation of alpha [0, inf) to [-1,1)'''
    return (1. - alpha ** 2.) ** 3. / (1. + alpha ** 3. + 4 * alpha ** 6.)


def LegVal(x: float, i: int) -> float:
    return float(np.polynomial.legendre.legval(x, [0] * i + [1]))


def LegVals(x: float, n: int) -> L[float]:
    '''len-n array of the Legendre polynomials evaluated at x'''
    return [LegVal(x, i) for i in range(n)]


# def LegProduct(s: float, a: float, A: np.ndarray = None) -> np.ndarray:
#     """
#     A - an MxN matrix with rows corresponding to s basis functions coefs
#         and columns corresponding to alpha basis function coefficients
#     Eq # 5 in mBEEF paper
#     """
#     out = [[ps*pa for pa in LegVals(t_alpha(a), 8)]
#            for ps in LegVals(t_s(s), 8)]
#     return np.array(out) if A is None else np.multiply(A, out)


def LegProduct(s: float, a: float, A: np.ndarray = None,
               trans: bool = False) -> np.ndarray:
    if A is None:
        coefs = np.ones((8, 8))
    else:
        coefs = A.reshape((8, 8))
    if trans:
        t_s, t_alpha = s, a
    else:
        p = s**2
        q = 0.804/(10./81.)
        t_s = 2.*p/(q+p) - 1.
        t_alpha = (1. - a ** 2.) ** 3. / (1. + a ** 3. + 4. * a ** 6.)
    ls, la = [LegVals(x, 8) for x in [t_s, t_alpha]]

    return np.multiply(coefs, np.outer(la, ls))
    # np.polynomial.legendre.legval2d(t_alpha, t_s, coefs)


# First derivatives

def dt_ds(s: float) -> float:
    '''First derivative of s transformation
    https://www.wolframalpha.com/input/
    ?i=first+derivative+of+2s%5E2%2F%28q%2Bs%29-1+with+respect+to+s
    '''
    if np.isinf(s):
        return 2.
    else:
        numer = 2 * s * (2*q + s)
        denom = (q + s)**2
        return numer/denom


def dt_da(alpha: float) -> float:
    n = 3*alpha*(8*alpha**4 + alpha**3 + alpha + 2)*(1-alpha**2)**2
    d = (4*alpha**6 + alpha**3 + 1)**2
    return n/d


def dL_dx(n: int, legvals: L[float]) -> Unary:
    '''First derivative of nth Legendre polynomial evalation.
    For efficiency, we rely on an array of precomputed values
    Li(x) at a specific value of x. Function is partially
    applied so that we can use it as a unary function.
    https://www.wolframalpha.com/input/
    ?i=first+derivative+of+LegendreP%5Bn%2C+x%5D
    '''
    def f(x: float) -> float:
        if n == 0:
            return 0
        elif abs(x) == 1:
            sign = 1. if x == 1 else (-1.)**(n+1)
            return sign * n*(n+1)/2
        elif x != 0:
            assert abs(x - x/abs(x)) > 1e-7
        numer = -(n+1)*(x*legvals[n] - legvals[n+1])
        denom = (x**2 - 1)
        return numer/denom
    return f


# Second derivatives

# def dt2_ds2(s: float) -> float:
#     '''Second derivative of s transformation
#     https://www.wolframalpha.com/input/
#     ?i=second+derivative+of+2s%5E2%2F%28q%2Bs%29-1+with+respect+to+s
#     '''
#     if np.isinf(s):
#         return 0.
#     else:
#         t1 = 4*s**2/(q+s)**3
#         t2 = 8*s/(q+s)**2
#         t3 = 4./(q+s)
#         return t1 - t2 + t3


# def dg2_ds2(a1: float = 4.9479) -> Unary:
#     def f(x: float) -> float:
#         if x == 0:
#             return 0
#         t1 = a1*float(np.exp(-a1*np.power(x, -0.5)))
#         t2 = 0.75 * x**3 - 0.25 * x**2.5
#         return t1*t2/x**5.5
#     return f


# def dt2_da2(x: float) -> float:
#     x2, x3, x4, x5, x6 = x**2, x**3, x**4, x**5, x**6
#     d = x6 + x3 + 1.
#     t1 = (6*(1-x2)**2 - 24*x2*(1-x2))/d
#     t2 = (12*x*(6.*x5 + 3*x2) * (1-x2)**2) / d**2
#     t31 = 2*(6*x5+3*x2)**2
#     t32 = 30*x4 + 6*x
#     t3 = (t31/d**3 - t32/d**2) * (1-x2)**3
#     return t1 - t2 - t3


def dL2_dx2(n: int, legvals: L[float]) -> Unary:
    '''Second derivative of nth Legendre polynomial evalation.
    For efficiency, we rely on an array of precomputed values
    Li(x) at a specific value of x. Function is partially
    applied so that we can use it as a unary function.

    https://www.wolframalpha.com/input/
    ?i=second+derivative+of+LegendreP%5Bn%2C+x%5D
    '''
    def f(x: float) -> float:
        if abs(x) == 1:
            sign = 1. if x == 1 else (-1.) ** n
            return sign * (n**2 + n - 2) * n * (n + 1)/8
        if x != 0:
            assert abs(x - x/abs(x)) > 1e-6
        x2 = x**2
        t1 = (n*(x2 - 1) + 2*x2)*legvals[n]
        t2 = 2*x*legvals[n+1]
        numer = t1 - t2
        denom = (x2 - 1)**2
        return (n+1) * numer / denom
    return f


# Third derivatives

# def dt3_ds3(s: float) -> float:
#     '''Third derivative of s transformation
#     https://www.wolframalpha.com/input/
#     ?i=third+derivative+of+2s%5E2%2F%28q%2Bs%29-1+with+respect+to+s
#     '''
#     if np.isinf(s):
#         return 0.
#     else:
#         t1 = 12*s**2/(q+s)**4
#         t2 = 24*s/(q+s)**3
#         t3 = 12/(q+s)**2
#         return t2 - (t1+t3)


# def dg3_dx3(a1: float = 4.9479) -> Unary:
#     def f(x: float) -> float:
#         if x == 0:
#             return 0
#         t1 = -a1/8*float(np.exp(-a1*np.power(x, -0.5)))
#         t2 = a1**2 - 9*a1*x**0.5 + 15*x
#         return t1*t2/(x**4.5)
#     return f


def dL3_dx3(n: int, legvals: L[float]) -> Unary:
    '''Third derivative of nth Legendre polynomial evalation.
    For efficiency, we rely on an array of precomputed values
    Li(x) at a specific value of x. Function is partially
    applied so that we can use it as a unary function.

    https://www.wolframalpha.com/input/
    ?i=third+derivative+of+LegendreP%5Bn%2C+x%5D
   '''
    def f(x: float) -> float:
        if abs(x) == 1:
            sgn = 1. if x == 1 else (-1.)**(n+1)
            return sgn * (n**4 + 2*n**3 - 7*n**2 - 8*n + 12) * n * (n + 1)/48
        assert abs(x - x/abs(x)) > 1e-5
        n2, n3 = n**2, n**3
        x2 = x**2
        t1 = (n+1) * (n2 * x2 - n2 + n*x2 - n + 6*x**2 + 2) * legvals[n+1]
        t2 = x * (n3*x2-n3+6*n2*x2-6*n2+11*n*x2-3*n+6*x2+2) * legvals[n]
        denom = (x2 - 1)**3
        return (t1-t2)/denom
    return f


# Inverses
def t_s_inv(trans_s: float) -> float:
    '''Inverse transformation: -1->1 to 0->inf.'''
    return -(q*(trans_s+1))**(0.5)/(trans_s-1)


def t_alpha_inv(a_: float) -> float:
    assert a_ <= 1 and a_ >= -.25
    res = opt.root_scalar(f=lambda x: t_alpha(x)-a_, x0=0, x1=10000)
    out = float(res.root)
    if np.isnan(out):
        pdb.set_trace()
    return out


def johannes_leg(n: int, x: float) -> np.ndarray:
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


if __name__ == '__main__':
    ss = np.linspace(0, 10, 1000)
    lvs = [LegVals(t_s(s), 7)for s in ss]
    l4s = [lv[5] for lv in lvs]
    import plotly
    import plotly.graph_objs as go
    data = [go.Scatter(x=ss, y=l4s, name='L5(x)')]
    f = go.Figure(data=data)
    plotly.offline.plot(f)
