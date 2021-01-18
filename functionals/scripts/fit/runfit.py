# External
from typing import Tuple as T, List as L


def runfit(data_: str, ce_: float, bm_: float, lc_: float, ref: str
           ) -> T[L[str], L[int], str, L[str], L[str], L[str],
                  L[str], L[str], L[str], L[int], str]:
    '''
    Arguments:
        data_ is serialized A*x<B matrices
        ce_, bm_, and lc_ are relative weights on the three datasets
            (lower = more weight)
        ref is the serialized functional that was used to generate the
            DFT data (i.e. the starting guess)

    Returns:
        xs_out: trajectory of functionals discovered by fitting
        opt_xs: indices into the above list - the subset that seem good
        (and other things)
    '''
    import json
    import functionals.fit.data as fdata
    import functionals.fit.functional as fx
    import functionals.scripts.fit.pareto as pareto
    import functionals.fit.math as math
    import plotly
    import plotly.graph_objs as go
    import numpy as np
    import scipy.optimize as opt
    from functionals.scripts.fit.legcountzeros import countzeros_in_s_at_alpha

    #############
    # Constants #
    #############
    n_cv = 0
    errtypes = ['ce', 'bm', 'lc', 'vol']
    num_maxiter = 300
    tol = 5e-6
    bounds = [(.5, 1.5)] + [(-.1, 0.1)] * 63
    rels = [False, True, True]
    #
    ms2fx = fx.FromMatrix.frompath('ms2')
    ms2 = ms2fx.x
    reffx = fx.FromMatrix(json.loads(ref))

    with open('/Users/ksb/scp_tmp/vi.json', 'r') as f:
        vi = json.load(f)
        viinv = np.linalg.inv(vi)

    leg_0 = np.array([
        -1., -1., -0.9491079123427593, -0.9324695142031523, -0.9061798459386646, -0.861136311594052, -0.7745966692414833, -0.7415311855993938, -0.6612093864662645, -0.5773502691896257, -0.5384693101056833, -0.40584515137739774, -0.3399810435848562, -0.238619186083197, -0.11, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

    midpoints = (leg_0[1:] + leg_0[:-1]) / 2.

    leg_0b = np.array([
        -0.25, -0.25, -0.238619186083197, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

    midpointsb = (leg_0b[1:] + leg_0b[:-1]) / 2.

    # apoints = [math.t_alpha(a) for a in np.linspace(0, 2, 10)]
    # spoints = [math.t_s(s) for s in np.linspace(0, 5.5, 55)]

    x_s = np.roll(range(len(midpoints)), len(midpoints) // 2 + 1
                  ) - float(int(len(midpoints) // 2))
    x_alpha = np.roll(range(len(midpointsb)), len(midpointsb) // 2 + 1
                      ) - float(int(len(midpointsb) // 2))
    m_x, m_y = np.meshgrid(x_s, x_alpha)
    l_s = np.log(x_s[1:len(x_s) // 2 + 1])
    l_alpha = np.log(x_alpha[1:len(x_alpha) // 2 + 1])
    l_x, l_y = np.meshgrid(l_s, l_alpha)

    # print("mpts based on initial guess, not ms2")
    # mpts = np.polynomial.legendre.leggrid2d(apoints, spoints,
    #                                         np.reshape(reffx.x, (8, 8))
    #                                         ).flatten()

    ################################################
    # Deserialize experimental and constraint data #
    ################################################

    data = fdata.Data.from_list(json.loads(data_))
    print('\nData len %d' % len(data))

    ce, bm, lc = refs = list(map(float, [ce_, bm_, lc_]))
    print('hyperparams are ', refs)

    As, Bs = data.xy()  # absolute CE, relative BM/LC
    Ax = np.vstack([a / w for a, w in zip(As, refs)])
    Bx = np.concatenate([b / w for b, w in zip(Bs, refs)])

    ################
    # Helper funcs #
    ################
    def data_penalty(x: np.ndarray) -> float:
        res = (Ax @ x) - Bx
        return float(np.sqrt(res @ res))

    def data_penalty_jac(x: np.ndarray) -> np.ndarray:
        res = (Ax @ x - Bx)
        num = res.T @ Ax
        denom = np.sqrt(res@res + 1)
        return num / denom

    def bound_penalty(x: np.ndarray) -> float:
        bcoeff = np.reshape(x, (8, 8))
        check = np.polynomial.legendre.leggrid2d(midpointsb, midpoints, bcoeff)
        negcheck = check[check < 0.]
        poscheck = check[check > 1.804]
        return float(np.sum(negcheck**2) + np.sum(poscheck**2))

    alphas = np.array([0, 0.5, 1, 1.5, 2, 100])
    ms2bmps = countzeros_in_s_at_alpha(ms2, alphas)
    bmpweight = np.array([10, 5, 3, 1, 0, 0, 0, 0])

    def count_bumps(x: np.ndarray) -> int:
        a0 = sum([bmpweight @ np.clip(np.array(bumpvec) - ms2bmp, 0, 10)
                  for ms2bmp, bumpvec
                  in zip(ms2bmps, countzeros_in_s_at_alpha(x, alphas))])
        return a0

    # def count_bumps(x: np.ndarray, debug: bool = False, lim: int = 0) -> int:
    #     '''Punish extra sign changes in 1st/2nd/3rd/4th derivative, relative
    #        to MS2'''
    #     bcoeff = np.reshape(x, (8, 8))
    #     check = np.polynomial.legendre.leggrid2d(apoints, spoints, bcoeff)
    #     penalty = []  # L[float]
    #     for deriv_bumps in [0, 1, 1, 2]:
    #         check = np.diff(check, axis=1)  # n'th derivative approximation
    #         bumps = np.abs(np.diff(np.sign(check))) / 2  # 1 for every bump
    #       penalty.append(sum([max(0, r.sum() - deriv_bumps) for r in bumps]))
    #         if deriv_bumps == 0:
    #             check = check[:, :60]
    #     return sum(penalty)

    # def curv_penalty(x: np.ndarray) -> float:
    #     '''Punish deviations from MS2'''
    #     bcoeff = np.reshape(x, (8, 8))
    #     check = np.polynomial.legendre.leggrid2d(apoints, spoints, bcoeff)
    #     res = check.flatten() - mpts
    #     return float(res @ res)

    def transform(x0: np.ndarray) -> np.ndarray:
        x1 = np.dot(viinv, x0)
        x1[:3] = [1., 6.696148885203612e+00, -2.351163592512346e-01]
        return np.dot(vi, x1)

    # def fitfun_(x0: np.ndarray, bumps: int, penalty: float) -> np.ndarray:
    #     def ff(z: np.ndarray) -> float:
    #         z = transform(z)
    #         bad = max(0., count_bumps(z) - bumps)
    #         bad += bound_penalty(z)
    #         bad += max(0., curv_penalty(z) - penalty)
    #         return 10000 * bad + data_penalty(z)
    #     sol = opt.minimize(ff, x0=x0,
    #                        method='nelder-mead', tol=tol,
    #                        options=dict(disp=True, maxiter=num_maxiter))
    #     out = transform(sol.x)
    #     return x0 if data_penalty(out) > data_penalty(x0) else out

    def fitfun(x0: np.ndarray, bumps: int, penalty: float) -> np.ndarray:

        cons = [
            opt.NonlinearConstraint(
                lb=-1, ub=bumps, keep_feasible=True,
                fun=lambda z: count_bumps(transform(z))),
            opt.NonlinearConstraint(
                lambda z: bound_penalty(transform(z)),
                lb=-1, ub=0.1, keep_feasible=True),
            # opt.NonlinearConstraint(
            #     lambda z: curv_penalty(transform(z)),
            #     lb=-1, ub=penalty, keep_feasible=True)
        ]

        sol = opt.minimize(
            fun=lambda z: data_penalty(transform(z)),
            jac=lambda z: data_penalty_jac(transform(z)),
            x0=x0, constraints=cons,
            method='trust-constr', tol=tol, bounds=bounds,
            options=dict(disp=2, maxiter=num_maxiter))
        out = transform(sol.x)
        return x0 if data_penalty(out) > data_penalty(x0) else out

    def show(x: np.ndarray, name: str = '') -> None:
        e = [data.mae(x, k, rel=r) for k, r in zip(errtypes, rels)]
        print('\n%s errs' % name, '  '.join(
            ['%s %.3f' % (k, e) for k, e in zip(errtypes, e)]))

    ##############
    # Initialize #
    ##############
    Lf = L[float]
    (xs, mcs, rmcs, mls, rmls, mvs, rmvs, mbs, rmbs, agg_err, curvs, bs) = (
        [], [], [], [], [], [], [], [], [], [], [], []
    )  # type: T[L[np.ndarray],Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf, L[int]]
    lists = [mcs, rmcs, mbs, rmbs, mls, rmls, mvs, rmvs]

    ###############################
    # Generate trajectory of fits #
    ###############################
    # ASSERT THAT THE INITIAL GUESS IS W/IN INITIAL CONSTRAINTS
    assert count_bumps(reffx.x) < 20
    # assert curv_penalty(reffx.x) < 1
    x = reffx.x
    for bump in [50, 100, 150]:
        # for curv_pen in [1, 2]:

        x = fitfun(x, bump, 0)
        xs.append(x)
        maes = [data.mae(x, k, r)
                for k in errtypes for r in [False, True]]
        print(bump, 0, maes)
        for lis, val in zip(lists, maes):
            lis.append(round(val, 3))
        agg_err.append(data_penalty(x))
        curvs.append(count_bumps(x))
        bs.append(count_bumps(x))

    #####################
    # VISUALIZE RESULTS #
    #####################
    datas, steps = [], []

    inds = [z[1] for z in pareto.pareto(
        [(cv, i) for i, cv in enumerate(curvs)], agg_err)[0]]
    inds = list(range(len(curvs)))  # comment out to only save pareto optimal

    for real_i, i in enumerate(sorted(inds)):
        p1, p2 = fx.FromMatrix(xs[i], name=str(i)).plot('black')
        datas.extend([p1, p2])
        tstr = 'Step %d curv %.2E err %.2f bump %d'
        title = tstr % (i, curvs[i], agg_err[i], bs[i])
        args = [dict(visible=[False] * 2 * len(inds)), {'title.text': title}]
        args[0]['visible'][2 * real_i:  # type: ignore
                           2 * real_i + 2] = [True, True]
        steps.append(dict(args=args, method='update'))

    curv_ = [curvs[i] for i in inds]
    agg_err_ = [agg_err[i] for i in inds]
    go.Figure(
        data=[go.Scatter(x=curv_, y=agg_err_, mode='markers',
                         hovertext=list(map(str, inds)))])
    fig2 = go.Figure(data=datas,
                     layout=dict(sliders=[dict(steps=steps)]))
    fig3 = go.Figure(
        data=[go.Scatter(x=curvs, y=agg_err, mode='lines+markers',
                         hovertext=list(map(str, range(len(xs)))))])

    for fig in [fig2, fig3]:
        plotly.offline.plot(fig)

    show(ms2, 'ms2')
    show(reffx.x, 'initial')
    show(x, 'final')

    opt_xs = list(range(len(xs)))  # all of them
    # list(map(int, input('\n\nWhich index is optimum?\n\n').split()))

    xs_out = [json.dumps(xs[i].tolist()) for i in opt_xs]
    bumps = [int(bs[i]) for i in opt_xs]

    stats = [[lst[i] for i in opt_xs] for lst in lists]

    #####################################
    # Cross validation - rewrite this? #
    #####################################
    cv = [dict(ce=[], bm=[], lc=[]),  # type: ignore
          dict(ce=[], bm=[], lc=[])]
    if n_cv > 0:
        for train, test in data.split(n_cv):
            # err = test.rmse if L2 else test.mae
            Atrn, Btrn = train.xy(rel=True)
            np.linalg.lstsq(Atrn, Btrn, rcond=None)[0]
            # xc = optim(Atrn, Btrn)
            # for i, xx in enumerate([xuc, xc]):
            #     for k in errtypes:
            #         cv[i][k].append(test.mae(xx, k))

    return (xs_out, opt_xs, json.dumps(cv),   # type: ignore
            *stats, None, None, bumps, '')

    # def l1(x: np.ndarray) -> float:
    #     '''scipy-cookbook.readthedocs.io/items/robust_regression.html'''
    #     res = (Ax @ x - Bx)
    #     return float(np.sum(np.sqrt(1 + np.square(res)) - 1))

    # def dl1(x: np.ndarray) -> np.ndarray:
    #     res = (Ax @ x - Bx)
    #     num = res.T @ Ax
    #     denom = np.sqrt(res@res + 1)
    #     return num / denom

    # def hess(x: np.ndarray) -> np.ndarray:
    #     return np.zeros((64, 64))

    # def nonlin(val: float) -> sp.optimize.NonlinearConstraint:
    #     return sp.optimize.NonlinearConstraint(
    #         fun=curv_penalty, lb=[0.], ub=[val], keep_feasible=False)

    # def optim(x_in: np.ndarray, cmax: float) -> np.ndarray:
    #     c = [constr.eq_constraint, constr.bound_con]  # , nonlin(cmax)]
    #     res = opt.minimize(
    #         x0=x, method='trust-constr', tol=TOL, constraints=c, fun=l1,
    # jac=dl1, hess=hess,
    #         options=dict(disp=True, maxiter=5000))
    #     return res.x

    # def fft_difference(x: np.ndarray) -> float:
    #     bcoeff = np.reshape(x, (8, 8))
    #     check = np.polynomial.legendre.leggrid2d(apoints, spoints, bcoeff)
    #     fft = np.log(np.abs(np.fft.fft2(check))[(m_x > 0) * (m_y > 0)])
    #     plane = np.linalg.lstsq(linreg2dX, fft, rcond=-1)
    #     SSres = np.sum(plane[1]**2)
    #     fftave = np.mean(fft)
    #     SStot = np.sum((fft - fftave)**2)
    #     Rsquared = 1. - SSres / SStot
    #     return float(1. - Rsquared)**2
