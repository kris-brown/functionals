from typing import List, Callable as C
import numpy as np
import json
import scipy.optimize as opt
from functionals.scripts.fit.count_bumps import count_bumps

# Constants
num_maxiter = 1000
tol = 5e-6
bounds = [(.5, 1.5)] + [(-.1, 0.1)] * 63

leg_0 = np.array([
    -1., -1., -0.9491079123427593, -0.9324695142031523, -0.9061798459386646, -0.861136311594052, -0.7745966692414833, -0.7415311855993938, -0.6612093864662645, -0.5773502691896257, -0.5384693101056833, -0.40584515137739774, -0.3399810435848562, -0.238619186083197, -0.11, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

midpoints = (leg_0[1:] + leg_0[:-1]) / 2.

leg_0b = np.array([
    -0.25, -0.25, -0.238619186083197, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

midpointsb = (leg_0b[1:] + leg_0b[:-1]) / 2.


with open('/Users/ksb/scp_tmp/vi.json', 'r') as f:
    vi = json.load(f)
    viinv = np.linalg.inv(vi)


def transform(x0: np.ndarray) -> np.ndarray:
    x1 = np.dot(viinv, x0)
    x1[:3] = [1., 6.696148885203612e+00, -2.351163592512346e-01]
    return np.dot(vi, x1)


def bound_penalty(x: np.ndarray) -> float:
    bcoeff = np.reshape(x, (8, 8))
    check = np.polynomial.legendre.leggrid2d(midpointsb, midpoints, bcoeff)
    negcheck = check[check < 0.]
    poscheck = check[check > 1.804]
    return float(np.sum(negcheck**2) + np.sum(poscheck**2))


def data_penalty(Ax: np.ndarray, Bx: np.ndarray) -> C[[np.ndarray], float]:
    def f(x: np.ndarray) -> float:
        res = (Ax @ x) - Bx
        return float(np.sqrt(res @ res))
    return f


def data_penalty_jac(Ax: np.ndarray, Bx: np.ndarray) -> C[[np.ndarray], float]:
    def f(x: np.ndarray) -> np.ndarray:
        res = (Ax @ x - Bx)
        num = res.T @ Ax
        denom = np.sqrt(res@res + 1)
        return num / denom
    return f


def fitfun(x0: np.ndarray, bumps: int,
           Ax: np.ndarray, Bx: np.ndarray,
           initial_bumps: List[List[int]]) -> np.ndarray:

    cons = [
        opt.NonlinearConstraint(
            lb=-1, ub=bumps, keep_feasible=True,
            fun=lambda z: count_bumps(transform(z), initial_bumps)),
        opt.NonlinearConstraint(
            lambda z: bound_penalty(transform(z)),
            lb=-1, ub=0.1, keep_feasible=True),
        # opt.NonlinearConstraint(
        #     lambda z: curv_penalty(transform(z)),
        #     lb=-1, ub=penalty, keep_feasible=True)
    ]
    dp = data_penalty(Ax, Bx)

    sol = opt.minimize(
        fun=lambda z: dp(transform(z)),
        jac=lambda z: data_penalty_jac(Ax, Bx)(transform(z)),
        x0=x0, constraints=cons,
        method='trust-constr', tol=tol, bounds=bounds,
        options=dict(disp=2, maxiter=num_maxiter))
    out = transform(sol.x)
    return x0 if dp(out) > dp(x0) else out
