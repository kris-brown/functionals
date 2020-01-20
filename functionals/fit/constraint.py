from typing import List, Tuple as T, Dict as D, Callable as C
import numpy as np
import scipy as sp
# Internal
import functionals.fit.math as math
import functionals.scripts.fit.h_norm_const as h_norm

#####################################


def h_norm_vec(_: float, __: float) -> np.ndarray:
    # a1 = 4.9479 msb = 1.
    return np.array(
        [h_norm.quad(h_norm.intkernel(i, j), 0., np.inf)[0]
         for j in range(8) for i in range(8)])


def L(i: int) -> math.Unary:
    return lambda ss: math.LegVal(math.t_s(ss), i)


def dL(i: int, legvals: List[float]) -> math.Unary:
    return math.chain_rule(math.dL_dx(i, legvals),
                           math.dt_ds, math.t_s)


def dL2(i: int, legvals: List[float]) -> math.Unary:
    return math.chain_rule2(math.dL2_dx2(i, legvals), math.dL_dx(i, legvals),
                            math.dt2_ds2, math.dt_ds, math.t_s)


def dL3(i: int, legvals: List[float]) -> math.Unary:
    return math.chain_rule3(
        math.dL3_dx3(i, legvals), math.dL2_dx2(i, legvals),
        math.dL_dx(i, legvals), math.dt3_ds3, math.dt2_ds2, math.dt_ds,
        math.t_s)


def deriv_constr(s: float, alpha: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the derivative of F(s,a) with respect to s.
    dF(s,a)/ds = Σ Σ xij * d(Li(t_s(s)))/ds * Lj(t_a(a)
    '''
    lv = math.LegVals(math.t_s(s), 9)  # precompute for efficiency
    # the values of Lj(alpha) pass through the derivative
    leg_alpha = math.LegVals(math.t_alpha(alpha), 8)
    # dg, g = math.dg_ds(a1), math.g_s(a1)

    vals = [math.chain_rule(dL(i, lv), math.dt_ds, math.t_s)(s)
            * leg_alpha[j]
            for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


def deriv2_constr(s: float, alpha: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the curvature of F(s,a) with respect to s.'''
    legvals = math.LegVals(math.t_s(s), 9)  # precompute for efficiency
    # dg2, dg, g = math.dg2_ds2(a1), math.dg_ds(a1), math.g_s(a1)

    # the values of Lj(alpha) pass through the derivative
    l_alpha = math.LegVals(math.t_alpha(alpha), 8)

    # d2F(s,a)/ds2 = Σ Σ xij * d2(Li(t_s(s)))/ds2 * Lj(t_a(a)
    vals = [math.chain_rule2(
        dL2(i, legvals), dL(i, legvals),
        math.dt2_ds2, math.dt_ds, math.t_s)(s)
        * l_alpha[j]
        for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


def deriv2_alpha(s: float, alpha: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the curvature of F(s,a) with respect to ALPHA.'''
    legvals = math.LegVals(math.t_alpha(alpha), 9)  # precompute for eff.
    # g = math.g_s(a1)(s)
    # the values of Lj(alpha) pass through the derivative
    legvals_s = math.LegVals(math.t_s(s), 8)

    # d2F(s,a)/da2 = Σ Σ xij *Li(t_s(s)) * d2(Lj(t_a(a))/da2
    vals = [legvals_s[i] * math.chain_rule2(
        math.dL2_dx2(i, legvals), math.dL_dx(i, legvals),
        math.dt2_da2, math.dt_da, math.t_alpha)(alpha)
        for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


def deriv3_constr(s: float, alpha: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the third derivative of F(s,a) with respect to s.'''
    lv = math.LegVals(math.t_s(s), 9)  # precompute for efficiency
    # dg3 = math.dg3_dx3(a1)
    # dg2, dg, g = math.dg2_ds2(a1), math.dg_ds(a1), math.g_s(a1)

    # the values of Lj(alpha) pass through the derivative
    legvals_alpha = math.LegVals(math.t_alpha(alpha), 8)
    # d3F(s,a)/ds3 = Σ Σ xij * d3(Li(t_s(s)))/ds3 * Lj(t_a(a)
    vals = [math.chain_rule3(
        dL3(i, lv), dL2(i, lv), dL(i, lv),
        math.dt3_ds3, math.dt2_ds2, math.dt_ds, math.t_s)(s)
        * legvals_alpha[j]
        for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


def point(s: float, a: float) -> np.ndarray:
    return math.LegProduct(s=s, a=a).flatten()


class Constraint(object):
    '''Sufficient information to *construct* a Ax=b matrix.
    Can either be POINT constraint or VECTOR constraint.'''
    GRID = 5

    def __init__(self, name: str, kind: str,
                 fun: C[[float, float], np.ndarray],
                 val: float, points: List[T[float, float]],
                 ) -> None:
        self.name = name
        self.kind = kind
        self.fun = fun
        self.points = points
        self.val = val
        assert kind in ['lt', 'gt']

    def __eq__(self, other: object) -> bool:
        return vars(self) == vars(other)

    def __len__(self) -> int:
        return len(self.points)

    def __str__(self) -> str:
        return '<%s>' % self.name

    __repr__ = __str__

    def err(self, x: np.ndarray) -> float:
        a, b = self.ab()
        return float(np.sum(a@x - b))

    def ab(self, points: List[T[float, float]] = None,
           cache_: D[T[float, float], List[float]] = None
           ) -> T[np.ndarray, np.ndarray]:
        A, vals = [], []  # type: T[List[List[float]], List[float]]
        cache = cache_ or {}

        neg = -1 if self.kind == 'gt' else 1
        for s, a in (self.points if points is None else points):
            # assert cache is None or (s,a) in cache,
            # (self.name, s, a, list(cache.keys()))
            pointval = self.fun(s, a) if cache is None else cache[(s, a)]
            A.append([x*neg for x in pointval])
            vals.append(self.val*neg)
        return np.array(A), np.array(vals)

    def cons(self, points: List[T[float, float]] = None,
             cache_: D[T[float, float], List[float]] = None
             ) -> T[np.ndarray, np.ndarray, np.ndarray]:
        a, ub = self.ab(points=points, cache_=cache_)
        lb = np.full(ub.shape, -np.inf)
        return a, lb, ub

    @staticmethod
    def lc(cons: List['Constraint']) -> sp.optimize.LinearConstraint:
        seen, reslt = set(),  []  # type: ignore
        for con in cons:
            a, ub = con.ab()
            for vec, u in zip(a, ub):
                key = tuple(vec.tolist())  # remove dups
                if key not in seen:
                    seen.add(key)
                    reslt.append((vec, -np.inf, u))
        alt, lblt, ublt = map(np.array, zip(*reslt))
        lclt = sp.optimize.LinearConstraint(alt, lblt, ublt)
        return lclt

    def test(self, pth: str, step: int = -1) -> None:
        '''Print constraint violations of a completed fit.'''
        from functionals.fit.fit import FitResult
        fr = FitResult.from_pth(pth)
        a, b = self.ab()
        x = fr.full.xs[step]
        res = a@x - b
        res = res[res > 0]  # if below zero, LT constraint is satisfied
        if len(res) == 0:
            res = [0]
        mean, median = np.mean(res), np.median(res)
        print('%s: %f %f' % (self.name, mean, median))


def mk_points(ts: float, ta: float, n_s: int, n_a: int,
              inv: bool = True) -> List[T[float, float]]:
    if inv:
        assert ts <= 1 and ta <= 1
        return [(math.t_s_inv(s), math.t_alpha_inv(a))
                for s in np.linspace(-1, ts, n_s)
                for a in np.linspace(-1, ta, n_a)]
    else:
        assert ts >= 0 and ta >= 0
        return [(s, a)
                for s in np.linspace(0, ts, n_s)
                for a in np.linspace(0, ta, n_a)]


leg_zeros = np.array([
    -1, -0.9491079123427593, -0.9324695142031523,
    -0.9061798459386646, -0.861136311594052, -0.7745966692414833,
    -0.7415311855993938, -0.6612093864662645, -0.5773502691896257,
    -0.5384693101056833, -0.40584515137739774, -0.3399810435848562,
    -0.238619186083197, -0.0, 0.23861918608319707, 0.3399810435848562,
    0.4058451513773972, 0.53846931010568, 0.57735026918962, 0.661209386466264,
    0.7415311855993948, 0.7745966692414835, 0.8611363115940525,
    0.9061798459386634, 0.9324695142031525, 0.9491079123427586, 1])

pos_vals = (leg_zeros[1:] + leg_zeros[:-1]) / 2   # MIDPOINTS

pos_points = [(math.t_s_inv(x), math.t_alpha_inv(y))
              for x in pos_vals for y in pos_vals]

ott = [1., 2, 3]
grid = [(s, a) for s in ott for a in ott]
shgrid = [(x-0.5, y-0.5) for x, y in grid]
szero = [(0., a) for a in ott]
azero = [(s, 0.) for s in ott]
dense = mk_points(4, 4, 20, 20, False)
constlist = [
    Constraint(name='lda_hi', fun=point, val=1.01, kind='lt', points=[(0, 1)]),
    Constraint(name='lda_lo', fun=point, val=.99, kind='gt', points=[(0, 1)]),
    Constraint(name='liebox', fun=point,  kind='lt', val=1.804,
               points=pos_points),
    Constraint(name='pos', fun=point, kind='gt', val=0.,
               points=pos_points),
    Constraint(name='zerocurv_hi', fun=deriv2_constr, points=[(0, 1)],
               val=20./81.+0.01, kind='lt'),
    Constraint(name='zerocurv_lo', fun=deriv2_constr, points=[(0, 1)],
               val=20./81.-0.01, kind='gt'),
    Constraint(name='curvneg', fun=deriv2_constr, kind='lt', val=1,
               points=dense),
    Constraint(name='curvpos', fun=deriv2_constr, kind='gt', val=-1,
               points=dense),
    Constraint(name='acurvneg', fun=deriv2_alpha, kind='lt', val=2,
               points=dense),  # mk_points(4, 1, 5, 5, False)),
    Constraint(name='acurvpos', fun=deriv2_alpha, kind='gt', val=-2,
               points=dense)]

#   Constraint(name='hnorm', fun=h_norm_vec, points=[(-1, -1)], kind='eq',
#              val=-.3125),
#   Constraint(name='scan11', fun=point, kind='lt', val=1.174,
#              ., points=azero),
#   Constraint(name='curv3', fun=deriv3_constr, kind='lt', val=-.1,
#              ., points=dense),
#   Constraint(name='decrease', fun=deriv_constr, kind='lt', val=0,
#              , points=[(s, a)
#                                 for s in np.linspace(2, 5, 10)
#                                 for a in np.linspace(0, 0.5, 5)]),
#   Constraint(name='increase', fun=deriv_constr, kind='gt', val=0,
#              , points=[(s, a)
#                                 for s in np.linspace(2, 5, 10)
#                                 for a in np.linspace(2, 3, 5)]),


# This constraint is trivial (vec is zero when applied at infinity)
# due to g(x) decay term
#   Constraint(name='decay', fun=deriv_constr, kind='lt', val=0.,
#             points=[(np.inf, a) for a in ott], default=0),

consts = {c.name: c for c in constlist}

if __name__ == '__main__':
    import pdb
    pdb.set_trace()
