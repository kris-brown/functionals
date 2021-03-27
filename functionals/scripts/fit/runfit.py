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
    import plotly
    import plotly.graph_objs as go
    import numpy as np
    from functionals.scripts.fit.count_bumps import count_bumps, all_bumps
    from functionals.scripts.fit.fitfun import fitfun, data_penalty
    ###############
    # Functionals #
    ###############
    ms2fx = fx.FromMatrix.frompath('ms2')
    mcaml = fx.FromMatrix.frompath('msurf')
    ms2 = ms2fx.x
    reffx = fx.FromMatrix(json.loads(ref))
    mcaml_mix = 0.0
    ###############
    # HYPERPARAMS #
    ###############
    errtypes = ['ce', 'bm', 'lc', 'vol']

    rels = [False, True, True]
    initial_guess = reffx.x * (1 - mcaml_mix) + mcaml.x * mcaml_mix

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

    initial_bumps = all_bumps(initial_guess)

    def show(x: np.ndarray, name: str = '') -> None:
        e = [data.mae(x, k, rel=r) for k, r in zip(errtypes, rels)]
        print('\n%s errs' % name, '  '.join(
            ['%s %.3f' % (k, e) for k, e in zip(errtypes, e)]))

    ##############
    # Initialize #
    ##############
    Lf = L[float]
    (xs, mcs, rmcs, mls, rmls, mvs, rmvs, mbs, rmbs, agg_err, bs) = (
        [], [], [], [], [], [], [], [], [], [], []
    )  # type: T[L[np.ndarray],Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf,Lf, L[int]]
    lists = [mcs, rmcs, mbs, rmbs, mls, rmls, mvs, rmvs]

    ###############################
    # Generate trajectory of fits #
    ###############################
    # ASSERT THAT THE INITIAL GUESS IS W/IN INITIAL CONSTRAINTS
    assert count_bumps(initial_guess, initial_bumps) < 10
    # assert curv_penalty(reffx.x) < 1
    x = initial_guess
    for bump in [5, 25]:
        # for curv_pen in [1, 2]:

        x = fitfun(x, bump, Ax, Bx, initial_bumps)
        xs.append(x)
        maes = [data.mae(x, k, r)
                for k in errtypes for r in [False, True]]
        print(bump, 0, maes)
        for lis, val in zip(lists, maes):
            lis.append(round(val, 3))
        agg_err.append(data_penalty(Ax, Bx)(x))
        bs.append(count_bumps(x, initial_bumps))

    #####################
    # VISUALIZE RESULTS #
    #####################
    datas, steps = [], []

    inds = [z[1] for z in pareto.pareto(
        [(cv, i) for i, cv in enumerate(bs)], agg_err)[0]]
    inds = list(range(len(bs)))  # comment out to only save pareto optimal

    for real_i, i in enumerate(sorted(inds)):
        p1, p2 = fx.FromMatrix(xs[i], name=str(i)).plot('black')
        datas.extend([p1, p2])
        tstr = 'Step %d err %.2f bump %d'
        title = tstr % (i, agg_err[i], bs[i])
        args = [dict(visible=[False] * 2 * len(inds)), {'title.text': title}]
        args[0]['visible'][2 * real_i:  # type: ignore
                           2 * real_i + 2] = [True, True]
        steps.append(dict(args=args, method='update'))

    curv_ = [bs[i] for i in inds]
    agg_err_ = [agg_err[i] for i in inds]
    go.Figure(
        data=[go.Scatter(x=curv_, y=agg_err_, mode='markers',
                         hovertext=list(map(str, inds)))])
    fig2 = go.Figure(data=datas,
                     layout=dict(sliders=[dict(steps=steps)]))
    fig3 = go.Figure(
        data=[go.Scatter(x=curv_, y=agg_err, mode='lines+markers',
                         hovertext=list(map(str, range(len(xs)))))])

    for fig in [fig2, fig3]:
        plotly.offline.plot(fig)

    show(ms2, 'ms2')
    show(initial_guess, 'initial')
    show(x, 'final')

    opt_xs = list(range(len(xs)))  # all of them
    # list(map(int, input('\n\nWhich index is optimum?\n\n').split()))

    xs_out = [json.dumps(xs[i].tolist()) for i in opt_xs]
    bumps = [int(bs[i]) for i in opt_xs]

    stats = [[lst[i] for i in opt_xs] for lst in lists]

    breakpoint()

    return (xs_out, opt_xs, "",   # type: ignore
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
