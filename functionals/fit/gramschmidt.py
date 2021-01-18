
from typing import Tuple
import numpy as np
from functionals.fit.constraint import point, curv
import functionals.scripts.fit.h_norm_const as h_norm


#################################################
aH = h_norm.hydro()
bH = -0.3125
aC = curv()
aC[0] += 0.0001
bC = 20./81.
aL = point(t_s=-1, t_a=0.)  # numerical issue if exactly 0
bL = 1
#################################################
AH = [-x/aH[0] for x in aH[1:]]
BH = bH/aH[0]


def h(p: np.ndarray, q: float) -> Tuple[np.ndarray, float]:
    '''Take 64 dim problem, eliminate first DOF'''
    return [p[0]*a + b for a, b in zip(AH, p[1:])], q-p[0]*BH


haC, hbC = h(aC, bC)

AC = [-x/haC[0] for x in haC[1:]]
BC = hbC/haC[0]


def c(p: np.ndarray, q: float) -> Tuple[np.ndarray, float]:
    '''Take 64 dim problem, eliminate first two DOF'''
    hp, hq = h(p, q)
    return [hp[0]*a+b for a, b in zip(AC, hp[1:])], hq - hp[0]*BC


caL, cbL = c(aL, bL)
AL = [-x/caL[0] for x in caL[1:]]
BL = cbL/caL[0]


def l(p: np.ndarray, q: float) -> Tuple[np.ndarray, float]:
    '''Take 64 dim problem, eliminate first three DOF'''
    cp, cq = c(p, q)
    return [cp[0]*a + b for a, b in zip(AL, cp[1:])], cq - cp[0]*BL


def compress(x: np.ndarray) -> np.ndarray:
    '''Take a 64 vector that may or may not satisfy the constraints and
    Compress it to 61 dimensions.'''
    A, b = transform(np.eye(64), x)
    y, resid, _, _ = np.linalg.lstsq(A, b, rcond=-1)
    if resid[0] > 0.01:
        print('WARNING: COMPRESSING TO 61 DOF WITHOUT SATISFYING CONSTRAINTS')

    return y


def transform(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    '''Apply h rowwise to the A matrix/B vector of a Ax=B'''
    tp, tq = map(np.array, zip(*[l(p, q) for p, q in zip(P, Q)]))
    return tp, tq


ms2 = [1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -
       2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06]


def recover(r: np.ndarray) -> np.ndarray:
    '''Take 61 dim vector and recover corresponding 64dim BEEF vector'''
    x3 = AL @ r + BL
    x2 = AC @ np.concatenate(([x3], r)) + BC
    x1 = AH @ np.concatenate(([x2, x3], r)) + BH
    return np.concatenate(([x1, x2, x3], r))


def consterr(x: np.ndarray) -> Tuple[float, float, float]:
    '''Report the error w/r/t constraints for a 64 length vector.'''
    constA = np.vstack([aL, aC, aH])
    constB = np.array([bL, bC, bH])
    a, b, c = (constA@x - constB).flatten()
    return a, b, c


def test() -> None:

    # Example: sum must equal 1
    A = np.ones((1, 64))
    B = np.ones(1)
    unconstrained = np.linalg.lstsq(A, B, rcond=-1)[0]

    unconst_err = consterr(unconstrained)

    tA, tB = transform(A, B)
    constrained_61 = np.linalg.lstsq(tA, tB, rcond=-1)[0]
    constrained_64 = recover(constrained_61)
    const_err = consterr(constrained_64)
    assert np.isclose(np.sum(unconstrained), 1)
    assert np.isclose(np.sum(constrained_64), 1)
    print('Err no constraint ', unconst_err, 'with constraint ', const_err)
