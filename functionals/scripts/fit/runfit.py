# External
from typing import Tuple as T, List as L


def runfit(data_: str, ce_: float, bm_: float, lc_: float, consts_: str,
           ) -> T[L[str], L[int], str, L[str], L[str], L[str],
                  L[str], L[str], L[str], str]:
    '''
    consts is serialized A*x<B matrices

    Returns res_vector, cv_data, mae_ce relmae_ce mae_bm relmae_bm
    mae_lat relmae_lat mae_vol relmae_vol mae_mag relmae_mag err
    '''
    import time
    import json
    import functionals.fit.data as fdata
    import functionals.fit.functional as fx
    import functionals.scripts.fit.pareto as pareto
    import functionals.fit.constraint as constr
    import plotly
    import plotly.graph_objs as go
    import pdb
    import numpy as np
    import scipy.optimize as opt
    import scipy as sp

    ms2 = x = np.array([1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -
                        2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06])

    ################
    # Helper funcs #
    ################

    # Curvature Penalty #
    s_vals, a_vals = constr.s_vals, constr.a_vals

    x_s = np.roll(list(range(len(s_vals))),
                  len(s_vals)//2+1) - float(int(len(s_vals)//2))
    x_a = np.roll(list(range(len(a_vals))),
                  len(a_vals)//2+1) - float(int(len(a_vals)//2))

    m_x, m_y = np.meshgrid(x_s, x_a)
    pos_x, pos_y = m_x[(m_x > 0) & (m_y > 0)], m_y[(m_x > 0) & (m_y > 0)]

    ms2grid = np.polynomial.legendre.leggrid2d(
        a_vals, s_vals, ms2.reshape((8, 8)))

    def curv_penalty(x: np.ndarray) -> np.ndarray:
        xx = x.reshape((8, 8))

        check_ = np.polynomial.legendre.leggrid2d(a_vals, s_vals, xx)
        check = check_ / (ms2grid**0.9)
        pos_check = np.abs(np.fft.fft2(check))[(m_x > 0) & (m_y > 0)]
        vander = np.column_stack((pos_x, pos_y, np.ones(len(pos_x))))
        sol = np.linalg.lstsq(vander, pos_check, rcond=-1)
        return sol[1]  # residual

    #############
    # Constants #
    #############
    n_cv = 0
    errtypes = ['ce', 'bm', 'lc', 'vol']
    TOL = 1e-10

    ##############
    # Initialize #
    ##############

    xs, mcs, rmcs, mls, rmls, mvs, rmvs, mbs, rmbs, agg_err, curvs = \
        [], [], [], [], [], [], [], [], [], [], []
    lists = [mcs, rmcs, mbs, rmbs, mls, rmls, mvs, rmvs]

    ################################################
    # Deserialize experimental and constraint data #
    ################################################

    data = fdata.Data.from_list(json.loads(data_))
    print('Data len %d' % len(data))
    # assert data.remove('ce', 'Cu')

    ce, bm, lc = refs = [0.05, 0.1, .1]  # list(map(float, [ce_, bm_, lc_]))

    As, Bs = data.xy(rel=True)
    Ax = np.vstack([a/w for a, w in zip(As, refs)])
    Bx = np.concatenate([b/w for b, w in zip(Bs, refs)])

    def l1(x: np.ndarray) -> float:
        '''scipy-cookbook.readthedocs.io/items/robust_regression.html'''
        res = (Ax @ x - Bx)
        return float(np.sum(np.sqrt(1+np.square(res))-1))

    def dl1(x: np.ndarray) -> np.ndarray:
        res = (Ax @ x - Bx)
        num = res.T @ Ax
        denom = np.sqrt(res@res + 1)
        return num / denom

    def hess(x: np.ndarray) -> np.ndarray:
        return np.zeros((64, 64))

    def nonlin(val: float) -> sp.optimize.NonlinearConstraint:
        return sp.optimize.NonlinearConstraint(
            fun=curv_penalty, lb=[0.], ub=[val], keep_feasible=False)  # 6.97 -> 6.9

    def optim(x_in: np.ndarray, cmax: float) -> np.ndarray:
        c = [constr.eq_constraint, constr.bound_con, nonlin(cmax)]
        res = opt.minimize(
            x0=x, method='trust-constr', tol=TOL, constraints=c,
            fun=l1, jac=dl1, hess=hess,
            options=dict(disp=False, maxiter=5000))
        return res.x

    print(curv_penalty(ms2))

    # xx = optim(ms2, np.inf)
    # print('ms2', [data.mae(ms2, k, rel=True) for k in errtypes])
    # print('overfit', [data.mae(xx, k, rel=True) for k in errtypes])
    # # go.Figure(data=fx.FromMatrix(xx, name='overfit').plot('black')).show()
    # import pdb; pdb.set_trace()
    # init_curv = curv_penalty(ms2)[0]
    # max_curv = curv_penalty(xx)[0]

    ###############################
    # Generate trajectory of fits #
    ###############################
    mcurvs = np.logspace(-6, 0, 30)
    mcurvs = np.linspace(6.1, 5.8, 3)
    ti = time.time()
    # mcurvs = np.logspace(-8, 1, 10)
    for i, curv in enumerate(mcurvs):
        print(i, round((time.time()-ti)/60))
        x = optim(x, curv)
        maes = [data.mae(x, k, r)
                for k in errtypes for r in [False, True]]
        for lis, val in zip(lists, maes):
            lis.append(round(val, 3))
        agg_err.append(rmcs[-1]/ce + rmbs[-1]/bm + rmls[-1]/lc)
        xs.append(x)
        curvs.append(curv_penalty(x)[0])
        print('\t thres: %.3E | curv: %.3E | err: %.2f'
              % (curv, curvs[-1], agg_err[-1]))

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
        args = [dict(visible=[False] * 2*len(inds)), {'title.text': title}]
        args[0]['visible'][2*real_i:2*real_i+2] = [True, True]  # type: ignore
        steps.append(dict(args=args, method='update'))

    curv_ = [curvs[i] for i in inds]
    agg_err_ = [agg_err[i] for i in inds]

    fig1 = go.Figure(
        data=[go.Scatter(x=curv_, y=agg_err_, mode='markers',
                         hovertext=list(map(str, inds)))])
    fig2 = go.Figure(data=datas,
                     layout=dict(sliders=[dict(steps=steps)]))
    fig3 = go.Figure(
        data=[go.Scatter(x=curvs, y=agg_err, mode='lines+markers',
                         hovertext=list(map(str, range(len(xs)))))])

    for fig in [fig1, fig2, fig3]:
        plotly.offline.plot(fig)

    fiterr = [data.mae(x, k, rel=True) for k in errtypes]
    finalerr = sum([x/y for x, y in zip(fiterr, refs)])
    print('\a\nfinal err %.2f' % finalerr, '  '.join(
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
            xuc = np.linalg.lstsq(Atrn, Btrn, rcond=None)[0]
            xc = optim(Atrn, Btrn)
            for i, xx in enumerate([xuc, xc]):
                for k in errtypes:
                    cv[i][k].append(test.mae(xx, k))

    return (xs_out, opt_xs, json.dumps(cv),   # type: ignore
            *stats, None, None, '')

    # lb, ub = map(list, zip(*bounds))

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

    # fun = fx.FromMatrix(res.x)
    # points = []  # type: L[L[float]]
    # success = True
    # cviol = round(l1const(res.x, A_con_, B_con_), 2)
    # for a in [0, 1, 2, 3]:
    #     points.append([])
    #     if not success:
    #         break
    #     for s in np.linspace(0, 4, 10):
    #         v = fun.apply(s, a)
    #         if v < -0.01 or v > 1.806:
    #             success = False
    #             print('\t failed bounds @(%.2f,%.2f) %.2f' % (s, a, v))
    #             break
    #         points[-1].append(v)
    # points = np.array(points)
    # if success:
    #     import pdb
    #     pdb.set_trace()

    # def l1(x: np.ndarray, A: np.ndarray, B: np.ndarray,
    #        A_lt: np.ndarray, B_lt: np.ndarray) -> float:
    #     '''Soft l1 loss
    #     https://scipy-cookbook.readthedocs.io/items/robust_regression.html'''
    #     res_expt = (A @ x - B) if np.size(B) else np.empty(0)
    #     res_lt = (A_lt @ x - B_lt) if np.size(B_lt) else np.empty(0)
    #     np.maximum(res_lt, 0, res_lt)  # fast relu (stackoverflow 32109319)
    #     res = np.concatenate([res_expt, res_lt])
    #     return float(np.sum(np.sqrt(1+np.square(res))-1))

    # def dl1(x: np.ndarray, A: np.ndarray, B: np.ndarray,
    #         A_lt: np.ndarray, B_lt: np.ndarray) -> np.ndarray:
    #     res_expt = (A @ x - B) if np.size(B) else np.empty(0)
    #     res_lt = (A_lt @ x - B_lt) if np.size(B_lt) else np.empty(0)
    #     np.maximum(res_lt, 0, res_lt)
    #     res = np.concatenate([res_expt, res_lt])
    #     num = res.T @ (np.vstack((A, A_lt)))
    #     denom = np.sqrt(res@res + 1)
    #     return num / denom

    # fun = fx.FromMatrix(gs.recover(x))
    # lines = np.array([[fun.apply(s, a)
    #                    for s in np.linspace(0, 4, 50)] for a in range(2)])
    # bends = np.max(
    #     np.sum(np.abs(np.diff(np.sign(np.diff(lines)))), axis=1))
    # if counter % 100 == 0:
    #     print(counter, bends)
    # if bends/2 > 2:
    #     fun.plot2d(True)
    #     pdb.set_trace()

    # ms2pow = np.ceil(np.log10(np.abs(raw_ms2)))
    # bounds = [(-a, a) for a in [10**(min(1, x+BOUND)) for x in ms2pow]]

    # initerr = l1expt(ms2, Ax61, Bx61)

    # # Convert these to 61 dimensional problems
    # # ------------------------------------------
    # Ax61, Bx61 = gs.transform(A_expt_, B_expt_)
    # Ac61, Bc61 = gs.transform(A_con_, B_con_)

    # CONSTRAINED OPTIMIZATION
    # --------------------------
    # overfit = np.linalg.lstsq(A_expt_, B_expt_, rcond=-1)[0]
    # overfiterr = [data.mae(overfit, k, rel=True) for k in errtypes]
    # tot_err = sum([x/y for x, y in zip(overfiterr, refs)])
    # print('\noverfit relative err %.2f' % tot_err, '  '.join(
    #     ['%s %.3f' % (k, e) for k, e in zip(errtypes, overfiterr)]))

    # overfit61 = np.linalg.lstsq(Ax61, Bx61, rcond=-1)[0]
    # overfit61err = [data.mae(gs.recover(overfit61), k, rel=True)
    #                 for k in errtypes]
    # tot_err61 = sum([x/y for x, y in zip(overfit61err, refs)])
    # print('\noverfit 61 relative err %.2f' % tot_err61, '  '.join(
    #     ['%s %.3f' % (k, e) for k, e in zip(errtypes, overfit61err)]))
    # overfit_cv = l1const(overfit61, Ac61, Bc61)

    # def loss(x: np.ndarray, As: L[np.ndarray], Bs: L[np.ndarray],
    #          weights: L[float],  # A_lt: np.ndarray, B_lt: np.ndarray,
    #          # x_weight: float = 0, init_curv: float = 0,
    #          debug: str = '') -> float:

    #     err_expt = sum(np.sum(np.abs(A @ x - B))/len(B)/w
    #                    for A, B, w in zip(As, Bs, weights))
    #     return
    # res_lt = (A_lt @ x - B_lt) if np.size(B_lt) else np.empty(0)
    # np.maximum(res_lt, 0, res_lt)
    # err_lt = np.sum(np.abs(res_lt))
    # err_curv = max(0, curv_penalty(x) - init_curv) * 1000
    # if debug:
    #     print("Debugging: "+debug+"\n\t with x_weight = %.2E" % x_weight)
    #     print(''.join(['\n\t%s %.2E' % (b, a) for a, b in zip(
    #         [err_expt, err_lt, err_curv],
    #         ['Err_expt', "Err_lt", "err_curv"])]))
    # return float(err_expt*x_weight + err_lt + err_curv)

    # ms2 = x = np.array([1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06])
