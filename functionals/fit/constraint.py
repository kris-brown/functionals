
from typing import List, Tuple as T, Dict as D, Callable as C
import numpy as np
import scipy as sp
import functools
# Internal
import functionals.fit.math as math
import functionals.scripts.fit.h_norm_const as h_norm

#####################################


@functools.lru_cache(maxsize=1000, typed=False)
def deriv_constr(t_s: float, t_a: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the derivative of F(ŝ,â) with respect to s.
    dF(ŝ,â)/dŝ = Σ Σ xij * dLi(ŝ)/dŝ * Lj(â)
    '''
    lv = math.LegVals(t_s, 9)  # precompute for efficiency
    # the values of Lj(alpha) pass through the derivative
    leg_alpha = math.LegVals(t_a, 8)

    vals = [math.dL_dx(i, lv)(t_s) * leg_alpha[j]
            for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


@functools.lru_cache(maxsize=1000, typed=False)
def deriv2_constr(t_s: float, t_a: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the curvature of F(ŝ,â) with respect to ŝ.'''
    legvals = math.LegVals(t_s, 9)  # precompute for efficiency

    # the values of Lj(alpha) pass through the derivative
    l_alpha = math.LegVals(t_a, 8)

    # d2F(ŝ,a)/ds2 = Σ Σ xij * d2(Li(ŝ))/dŝ2 * Lj(t_a(a)
    vals = [math.dL2_dx2(i, legvals)(t_s) * l_alpha[j]
            for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


@functools.lru_cache(maxsize=1000, typed=False)
def deriv2_alpha(t_s: float, t_a: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the curvature of F(s,a) with respect to ALPHA.'''
    legvals = math.LegVals(t_a, 9)  # precompute for eff.
    # g = math.g_s(a1)(s)
    # the values of Lj(alpha) pass through the derivative
    legvals_s = math.LegVals(t_s, 8)

    # d2F(s,a)/da2 = Σ Σ xij *Li(t_s(s)) * d2(Lj(t_a(a))/da2
    vals = [legvals_s[i] * math.dL2_dx2(j, legvals)(t_a)
            for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


def deriv3_constr(t_s: float, t_a: float) -> np.ndarray:
    '''Return a vector which, when dotted with a BEEF vector,
    gives the third derivative of F(s,a) with respect to s.'''
    lv = math.LegVals(t_s, 9)  # precompute for efficiency

    # the values of Lj(alpha) pass through the derivative
    legvals_alpha = math.LegVals(t_a, 8)
    # d3F(s,a)/ds3 = Σ Σ xij * d3(Li(t_s(s)))/ds3 * Lj(t_a(a)
    vals = [math.dL3_dx3(i, lv)(t_s) * legvals_alpha[j]
            for j in range(8) for i in range(8)]
    return np.array(vals)  # flattened 64-element array


@functools.lru_cache(maxsize=1000, typed=False)
def point(t_s: float, t_a: float) -> np.ndarray:
    return math.LegProduct(s=t_s, a=t_a, trans=True).flatten()


class Constraint(object):
    '''Sufficient information to *construct* a Ax=b matrix.
    Can either be POINT constraint or VECTOR constraint.'''

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
        return float(np.sum(np.maximum(a@x - b, 0)))

    def ab(self) -> T[np.ndarray, np.ndarray]:
        A, vals = [], []  # type: T[List[List[float]], List[float]]
        # cache = cache_ or {}

        neg = -1 if self.kind == 'gt' else 1
        for s, a in self.points:
            pointval = self.fun(s, a)  # if cache is None else cache[(s, a)]
            assert not np.isnan(pointval).any()
            assert not np.isinf(pointval).any()
            A.append([x*neg for x in pointval])
            vals.append(self.val*neg)
        return np.array(A), np.array(vals)

    @staticmethod
    def AB(cons: D['str', int]) -> T[np.ndarray, np.ndarray]:
        '''Produce vstacked A,b of many constraints.'''
        A, B = np.empty((0, 64)), np.empty(0)
        for c in constlist:
            weight_ = cons.get(c.name, 0)
            weight = 10**weight_
            a, b = map(np.array, c.ab())
            A = np.vstack((A, a*weight))
            B = np.concatenate((B, b*weight))
        return A, B


leg_zeros = np.array(np.array([-1., -1., -0.9491079123427593, -0.9324695142031523, -0.9061798459386646, -0.861136311594052, -0.7745966692414833, -0.7415311855993938, -0.6612093864662645, -0.5773502691896257, -0.5384693101056833, -0.40584515137739774, -0.3399810435848562, -0.238619186083197, -0.11, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])
                     )

s_vals = (leg_zeros[1:] + leg_zeros[:-1]) / 2   # MIDPOINTS
a_vals = [-.25]+[a for a in s_vals if a > -0.25]
pos_points = [(x, y) for x in s_vals for y in a_vals]
azero = [(s, 1.) for s in s_vals]
sinf = [(1., float(a)) for a in a_vals]

dense = [(math.t_s(s), math.t_alpha(a))
         for s in np.concatenate((np.linspace(0, 5, 10), [np.inf]))
         for a in np.linspace(0, 10, 10)]

constlist = [
    Constraint(name='liebox', fun=point,  kind='lt', val=1.804,
               points=pos_points),
    Constraint(name='pos', fun=point, kind='gt', val=0.,
               points=pos_points),
    # CANNOT ADD THESE UNTIL WE INITIALIZE WITH BEEF-SCAN
    # Constraint(name='scan11', fun=point, kind='lt', val=1.174, points=azero),
    # Constraint(name='decay', fun=point, kind='lt', val=0.1, points=sinf),
]

# Constraint(name='curvneg', fun=deriv2_constr, kind='lt', val=0.3*1.1,
# Constraint(name='curvpos', fun=deriv2_constr, kind='gt', val=-0.3*1.1,
# Constraint(name='acurvneg', fun=deriv2_alpha, kind='lt', val=.002*1.1,
# Constraint(name='acurvpos', fun=deriv2_alpha, kind='gt', val=-.002*1.1,
# Constraint(name='curv3pos', fun=deriv3_constr, kind='lt', val=.3*1.1,
# Constraint(name='curv3neg', fun=deriv3_constr, kind='gt', val=-.3*1.1,
#   Constraint(name='increase', fun=deriv_constr, kind='gt', val=0,
#              , points=[(s, a) for s in np.linspace(2, 5, 10)
#                               for a in np.linspace(2, 3, 5)]),


def curv() -> np.ndarray:
    # weights curvature about s=0,a=1
    talpha1 = 0.
    q = 0.804/(10./81.)
    leg8talpha1 = math.johannes_leg(8, talpha1)
    n8 = np.linspace(0., 7., 8)
    curvs0 = -(-1.)**(n8) * n8*(n8+1.)/2. * 4./q
    return np.outer(leg8talpha1, curvs0).flatten()


consts = {c.name: c for c in constlist}

eq_A = [point(-1., 0.), curv(), h_norm.hydro()]
eq_B = [1, 20./81., -.3125]
eq_constraint = sp.optimize.LinearConstraint(eq_A, eq_B, eq_B,
                                             keep_feasible=True)

quasi_eq = [sp.optimize.LinearConstraint(A, b-0.05, b + 0.05,
                                         keep_feasible=True)
            for A, b in zip(eq_A, eq_B)]

bound_con = sp.optimize.LinearConstraint(
    A=[point(s, a) for s, a in pos_points],
    lb=np.zeros(len(pos_points)),
    ub=np.ones(len(pos_points))*1.804,
    keep_feasible=True)
bound_con_dense = sp.optimize.LinearConstraint(
    A=[point(s, a) for s, a in dense],
    lb=np.zeros(len(dense)),
    ub=np.ones(len(dense))*1.804,
    keep_feasible=True)

ms2 = [1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -
       2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06]

if __name__ == '__main__':
    pass
    # import functionals.fit.functional as fx

    # s = 3
    # ts = math.t_s(s)
    # lv = math.LegVals(ts, 9)
    # print(math.dL_dx(0, lv)(ts))
    # ss = np.linspace(0, 5, 500)
    # x = np.array([0]*12+[1]+[0]*51)
    # fun = fx.FromMatrix(x)
    # ys = [fun.apply(s, 1.5) for s in ss]
    # dys = [deriv_constr(s, 1.5)@x for s in ss]
    # data = [go.Scatter(x=ss, y=ys, name='L5(x)', mode='markers'),
    #         go.Scatter(x=ss, y=dys, name='dL5(x)/dx', mode='markers')]
    # f = go.Figure(data=data)
    # plotly.offline.plot(f)

    # fx.plots([fx.MS2, fx.FromMatrix(ms2)], True)
