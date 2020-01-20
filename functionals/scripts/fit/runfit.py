# External
from typing import Tuple as T


def runfit(data_: str, ce_: float, bm_: float, lc_: float, consts_: str,
           condata: str) -> T[str, str, float, float, float]:
    '''
    consts is a dict from const name to list of points
        (if different from default)
    condata is map from name to:
        (ABdata, kind, dfaultval, default_points)
    '''

    import json
    import ast
    # import functionals.fit.functional as fx
    import functionals.fit.data as fdata
    import functionals.fit.constraint as constr
    import numpy as np
    import scipy.optimize as opt

    # Constants
    L2 = False
    do_cv = False
    n_cv = 5
    pbesol = np.array([1.402, 0.402]+[0]*62)

    # Cost functions

    def l2(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> float:
        res = A @ x - B
        return float(res @ res)

    def dl2(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        res = A @ x - B
        return 2 * res.T @ A

    def l1(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> float:
        '''Soft l1 loss
        https://scipy-cookbook.readthedocs.io/items/robust_regression.html'''
        res = A @ x - B
        res2 = res @ res
        return float(np.sum(np.sqrt(1+res2)-1))

    def dl1(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        res = A @ x - B
        num = res.T @ A
        denom = np.sqrt(res@res + 1)
        return num / denom

    fun = dict(fun=l2, jac=dl2) if L2 else dict(fun=l1, jac=dl1)

    # Bounds
    beefpow = [1., 0., 0., -1., -2., -1., -2., -2.,
               0., 0., -1., -2., -1., -2., -2., -2.,
               -1., -1., -1., -1., -1., -1., -6., -6.,
               -1., -2., -2., -2., -2., -2., -7., -7.,
               -2., -2., -2., -2., -6., -6., -7., -7.,
               -3., -3., -3., -3., -7., -7., -7., -7.,
               -4., -4., -7., -7., -7., -7., -8., -8.,
               -5., -5., -8., -8., -8., -8., -8., -8.]
    bounds = [(-a, a) for a in [10**(min(1, x+10)) for x in beefpow]]
    lb, ub = map(list, zip(*bounds))

    # Deserialize experimental and constraint data

    data = fdata.Data.from_list(json.loads(data_))
    err = data.rmse if L2 else data.mae
    ce, bm, lc = map(float, [ce_, bm_, lc_])
    bm = lc = 1000000000

    A_expt, B_expt = [z*100 for z in data.xy(ce, bm, lc)]  # change scale

    consts = json.loads(consts_)
    cache = {k: v[0] for k, v in ast.literal_eval(condata).items()}

    # ASSEMBLE A/B FROM ABDATA and POINTS
    ckwargs = [dict(points=consts.get(c.name), cache_=cache[c.name])
               for c in constr.constlist]

    # MAKE CONSTRAINTS

    As, Ls, Us = zip(*filter(lambda x: x[0].size, [
        c.cons(**k) for c, k in zip(constr.constlist, ckwargs)]))

    lclt = opt.LinearConstraint(np.vstack(As), np.concatenate(Ls),
                                np.concatenate(Us))

    constraints = []  # [lclt]

    # CONSTRAINED OPTIMIZATION

    overfit = np.linalg.lstsq(A_expt, B_expt, rcond=-1)[0]
    ms2 = [1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -
           2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06]

    opt_args = dict(
        method='slsqp', x0=ms2, bounds=bounds, tol=1e-10,
        constraints=constraints,  **fun,
        options=dict(disp=True, maxiter=5000))

    x = opt.minimize(args=(A_expt, B_expt), **opt_args).x

    # Error analysis
    errtypes = ['ce', 'bm', 'lc']
    mc, mb, ml = [err(x, k) for k in errtypes]

    print(mc, mb, ml)
    import pdb
    pdb.set_trace()

    # Cross validation
    cv = [dict(ce=[], bm=[], lc=[]),  # type: ignore
          dict(ce=[], bm=[], lc=[])]
    if do_cv:
        for train, test in data.split(n_cv):
            err = test.rmse if L2 else test.mae
            Atrn, Btrn = train.xy(ce, bm, lc)
            xuc = np.linalg.lstsq(Atrn, Btrn, rcond=None)[0]
            xc = opt.minimize(args=(Atrn, Btrn), **opt_args).x
            for i, xx in enumerate([xuc, xc]):
                for k in errtypes:
                    cv[i][k].append(err(xx, k))

    return json.dumps(x.tolist()), json.dumps(cv), mc, mb, ml

    # Test the unconstrained lst sq solution
    # xx = np.linalg.lstsq(A_expt, B_expt, rcond=None)[0]
    # mc, mb, ml = [data.mae(xx, k) for k in ['ce', 'bm', 'lc']]
    # return json.dumps(xx.tolist()), '[]', mc, mb, ml

    # ABeqs = [c.cons_eq(**k) for c, k in zip(constr.constlist, ckwargs)]
    # A_eqs, B_eqs = zip(*filter(lambda x: x[0].size, ABeqs))

    # A_const = np.vstack(A_eqs)
    # B_const = np.concatenate(B_eqs)

    # x0 = opt.least_squares(
    #     cost, fx.BEEF().x, jac=dcost, bounds=(lb, ub),
    #     kwargs={'A': A_const, 'B': B_const}).x

    # seen, reseq, reslt = set(), [], []  # type: ignore
    # for a, ub, kind, _ in cons:
    #     for vec, u in zip(a, ub):
    #         key = tuple(vec)  # remove dups
    #         if key not in seen:
    #             seen.add(key)
    #             if kind == 'eq':
    #                 reseq.append((vec, u, u))
    #             else:
    #                 reslt.append((vec, -np.inf, u))
    # aeq, lbeq, ubeq = map(np.array, zip(*reseq))
    # alt, lblt, ublt = map(np.array, zip(*reslt))
    # clt = opt.LinearConstraint(alt, lblt, ublt)

    def cost(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return A @ x - B

    def dcost(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return A

#   beef = np.array([1.18029330e+00, 8.53027860e-03, -1.02312143e-01, 6.85757490e-02, -6.61294786e-03, -2.84176163e-02, 5.54283363e-03, 3.95434277e-03, -1.98479086e-03, 1.00339208e-01, -4.34643460e-02, -1.82177954e-02, 1.62638575e-02, -8.84148272e-03, -9.57417512e-03, 9.40675747e-03, 6.37590839e-03, -8.79090772e-03, -1.50103636e-02, 2.80678872e-02, -1.82911291e-02, -1.88495102e-02, 1.69805915e-07, -2.76524680e-07, 1.44642135e-03, -3.03347141e-03, 2.93253041e-03, -8.45508103e-03, 6.31891628e-03, -8.96771404e-03, -2.65114646e-08, 5.05920757e-08,
#                    6.65511484e-04, 1.19130546e-03, 1.82906057e-03, 3.39308972e-03, -7.90811707e-08, 1.62238741e-07, -4.16393106e-08, 5.54588743e-08, -1.16063796e-04, 8.22139896e-04, -3.51041030e-04, 8.96739466e-04, 2.09603871e-08, -3.76702959e-08, 2.36391411e-08, -3.38128188e-08, -5.54173599e-06, -5.14204676e-05, 6.68980219e-09, -2.16860568e-08, 9.12223751e-09, -1.38472194e-08, 6.94482484e-09, -7.74224962e-09, 7.36062570e-07, -9.40351563e-06, -2.23014657e-09, 6.74910119e-09, -4.93824365e-09, 8.50272392e-09, -6.91592964e-09, 8.88525527e-09])
