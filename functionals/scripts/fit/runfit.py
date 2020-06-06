# External
from typing import Tuple as T, List as L


def runfit(data_: str, ce_: float, bm_: float, lc_: float, consts_: str
           ) -> T[L[str], L[int], str, L[str], L[str], L[str],
                  L[str], L[str], L[str], str]:
    '''
    consts is serialized A*x<B matrices

    Returns res_vector, cv_data, mae_ce relmae_ce mae_bm relmae_bm
    mae_lat relmae_lat mae_vol relmae_vol mae_mag relmae_mag err
    '''
    import json
    import functionals.fit.data as fdata
    import functionals.fit.functional as fx
    import functionals.scripts.fit.pareto as pareto
    import plotly
    import plotly.graph_objs as go
    import pdb
    import numpy as np
    import scipy.optimize as opt

    ms2 = x = np.array([
        1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06])
    #############
    # Constants #
    #############
    n_cv = 0
    errtypes = ['ce', 'bm', 'lc', 'vol']

    boundpenalty = 10
    num_maxiter = 100000
    tol = 1e-7

    with open('/Users/ksb/scp_tmp/vi.json', 'r') as f:
        vi = json.load(f)
    viinv = np.linalg.inv(vi)

    leg_0 = np.array([
        -1., -1., -0.9491079123427593, -0.9324695142031523, -0.9061798459386646, -0.861136311594052, -0.7745966692414833, -0.7415311855993938, -0.6612093864662645, -0.5773502691896257, -0.5384693101056833, -0.40584515137739774, -0.3399810435848562, -0.238619186083197, -0.11, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

    midpoints = (leg_0[1:] + leg_0[:-1]) / 2.

    leg_0b = np.array([
        -0.25, -0.25, -0.238619186083197, -0.0, 0.23861918608319707, 0.3399810435848562, 0.4058451513773972, 0.538469310105683, 0.577350269189626, 0.6612093864662646, 0.7415311855993948, 0.7745966692414835, 0.8611363115940525, 0.9061798459386636, 0.9324695142031525, 0.9491079123427586, 1.0, 1.0])

    midpointsb = (leg_0b[1:] + leg_0b[:-1]) / 2.

    x_s = np.roll(range(len(midpoints)), len(midpoints) // 2 + 1
                  ) - float(int(len(midpoints) // 2))

    x_alpha = np.roll(range(len(midpointsb)), len(midpointsb) // 2 + 1
                      ) - float(int(len(midpointsb) // 2))

    m_x, m_y = np.meshgrid(x_s, x_alpha)

    l_s = np.log(x_s[1:len(x_s) // 2 + 1])
    l_alpha = np.log(x_alpha[1:len(x_alpha) // 2 + 1])
    l_x, l_y = np.meshgrid(l_s, l_alpha)

    llx = l_x.flatten()
    lly = l_y.flatten()

    linreg2dX_ = np.array([llx, lly]).T
    linreg2dX = np.c_[linreg2dX_, np.ones(linreg2dX_.shape[0])]

    ################################################
    # Deserialize experimental and constraint data #
    ################################################

    data = fdata.Data.from_list(json.loads(data_))
    print('\nData len %d' % len(data))

    ce, bm, lc = refs = [1, 8, 4]  # list(map(float, [ce_, bm_, lc_]))

    As, Bs = data.xy(rel=True)
    Ax = np.vstack([a / w for a, w in zip(As, refs)])
    Bx = np.concatenate([b / w for b, w in zip(Bs, refs)])

    ################
    # Helper funcs #
    ################
    def data_penalty(x: np.ndarray) -> float:
        '''scipy-cookbook.readthedocs.io/items/robust_regression.html'''
        res = (Ax @ x) - Bx
        return float(np.sqrt(res @ res))  # float(np.sum(np.sqrt(1 + np.square(res)) - 1))  # L1 norm

    def bound_penalty(x: np.ndarray) -> float:
        bcoeff = np.reshape(x, (8, 8))
        check = np.polynomial.legendre.leggrid2d(midpointsb, midpoints, bcoeff)
        negcheck = check[check < 0.]
        poscheck = check[check > 1.804]
        return float(np.sum(negcheck**2) + np.sum(poscheck**2))

    def curv_penalty(x: np.ndarray) -> float:
        # Curvature Penalty #
        bcoeff = np.reshape(x, (8, 8))
        check = np.polynomial.legendre.leggrid2d(midpointsb, midpoints, bcoeff)
        fft = np.log(np.abs(np.fft.fft2(check))[(m_x > 0) * (m_y > 0)])
        plane = np.linalg.lstsq(linreg2dX, fft, rcond=-1)

        SSres = np.sum(plane[1]**2)
        fftave = np.mean(fft)
        SStot = np.sum((fft - fftave)**2)
        Rsquared = 1. - SSres / SStot
        return float(1. - Rsquared)**2

    def transform(x0: np.ndarray) -> np.ndarray:
        x1 = np.dot(viinv, x0)
        x1[:3] = [1., 6.696148885203612e+00, -2.351163592512346e-01]
        return np.dot(vi, x1)

    def fitfun_(x0: np.ndarray, fftpenalty: float) -> float:
        x = transform(x0)

        # residual with data
        datares = data_penalty(x)
        # inequality constraints
        bcheck = boundpenalty * bound_penalty(x)
        # Curvature constraints
        fftres = fftpenalty * curv_penalty(x)

        return float(datares + fftres + bcheck)

    def fitfun(x0: np.ndarray, fftpenalty: float) -> np.ndarray:
        sol = opt.minimize(fun=fitfun_, x0=x, args=(fftpenalty,),
                           method='nelder-mead', tol=tol,
                           options=dict(disp=True, maxiter=num_maxiter))
        return transform(sol.x)
    ##############
    # Initialize #
    ##############
    Lf = L[float]
    (xs, mcs, rmcs, mls, rmls, mvs, rmvs, mbs, rmbs, agg_err, curvs) = (
        [], [], [], [], [], [], [], [], [], [], []
    )  # type: T[L[np.ndarray],Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf]
    lists = [mcs, rmcs, mbs, rmbs, mls, rmls, mvs, rmvs]

    ###############################
    # Generate trajectory of fits #
    ###############################

    for curv_pen in np.logspace(0, 5, 5):
        x = fitfun(ms2, curv_pen)
        xs.append(x)

        maes = [data.mae(x, k, r)
                for k in errtypes for r in [False, True]]
        for lis, val in zip(lists, maes):
            lis.append(round(val, 3))
        agg_err.append(data_penalty(x))
        curvs.append(curv_penalty(x))

    #####################
    # VISUALIZE RESULTS #
    #####################
    datas, steps = [], []
    inds = [z[1] for z in pareto.pareto(
        [(cv, i) for i, cv in enumerate(curvs)], agg_err)[0]]

    for real_i, i in enumerate(sorted(inds)):
        p1, p2 = fx.FromMatrix(xs[i], name=str(i)).plot('black')
        datas.extend([p1, p2])
        title = 'Step %d curv %.2E err %.2f' % (i, curvs[i], agg_err[i])
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

    initerr = [data.mae(ms2, k, rel=True) for k in errtypes]
    fiterr = [data.mae(x, k, rel=True) for k in errtypes]
    print('\ninitial errs', '  '.join(
        ['%s %.3f' % (k, e) for k, e in zip(errtypes, initerr)]))
    print('\a\nfinal errs', '  '.join(
        ['%s %.3f' % (k, e) for k, e in zip(errtypes, fiterr)]))

    opt_xs = list(map(int, input('\n\nWhich index is optimum?\n\n').split()))

    xs_out = [json.dumps(xs[i].tolist()) for i in opt_xs]
    stats = [[lst[i] for i in opt_xs] for lst in lists]

    pdb.set_trace()

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
            *stats, None, None, '')

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
