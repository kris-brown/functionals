# External modules
from typing import Tuple as T, Any, Optional, Dict, List
import ast
import collections
import csv
import itertools
import os
import copy
import json
import sys
import numpy as np
import plotly
import plotly.graph_objs as go
import ase.units as units
import ase.io as aseio
import ase
import subprocess as subp
# Internal Modules
import dbgen
import functionals.fit.functional as fx
import functionals.fit.data as fdata
from functionals.CLI.submit import matdata
from functionals.scripts.io.parse_job import parse_job
from functionals.scripts.fit.pareto import pareto

"""
A CLI interface for common queries to the DB
"""


#############################################################################
datapth = '/' + os.path.join(*__file__.split('/')[:-3], 'data/')
volcsv = datapth + 'opt_vol.csv'
db = datapth + 'functionals.json'
default = dbgen.ConnectInfo.from_file(db)
ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3


def safe_sub(x: Optional[float], y: Optional[float]) -> Optional[float]:
    if None in [x, y]:
        return None
    else:
        assert x and y
        return x - y


def mbfloat(x: str) -> Optional[float]:
    return float(x) if x else None


def safeAvg(xs: List[Optional[float]]) -> Optional[float]:
    notnull = list(filter(None, xs))
    # if any([isinstance(x, str) for x in xs]):
    #     print(xs); ; breakpoint();
    return sum(notnull) / len(notnull) if notnull else None

###############
# DFT RELATED #
###############


def viz_plts(name: str, calc: str, v_: str, e_: str,
             bulkmod_: str, vol_: str, eng_: str, abc: str,
             expt_vol: float, expt_bm: float,
             bm_ls: float) -> T[str, Any, Any, Any]:
    v, e = map(json.loads, [v_ or '[]', e_ or '[]'])
    bulkmod, vol, eng = [float(x) if x else 0. for x in [bulkmod_, vol_, eng_]]
    m, m2, m3 = [dict(color='rgb(%d, %d, %d)' % (r, g, b), size=s)
                 for r, g, b, s
                 in [(231, 99, 250, 20), (99, 231, 250, 5), (199, 231, 50, 5)]]

    title = '''{} {}<br><b>EOS</b>: vol={:.1f} Å³,<br><b>Expt</b>: vol = {:.1f} Å³, bm = {:.1f}
                GPa<br><b>Stencil</b> bm: {:.1f} GPa,
                <b>QuadFit</b> bm: {:.1f} GPa'''
    title = title.format(name, calc, vol, expt_vol or 0,
                         expt_bm or 0, bulkmod, bm_ls or 0)
    if abc:
        a, b, c = json.loads(abc)
        quad_v = np.linspace(min(v), max(v))
        quad_y = [a * x**2 + b * x + c for x in quad_v]
        bulkmod_a = bulkmod / (2 * vol * ev_a3_to_gpa)
        bulkmod_y = [bulkmod_a * (x - vol)**2 + eng for x in quad_v]
    else:
        quad_v = quad_y = bulkmod_y = []
    return (title,
            go.Scatter(x=v, y=e, mode='markers', marker=m, name='opt'),
            go.Scatter(x=quad_v, y=quad_y, mode='markers', marker=m2,
                       name='quad'),
            go.Scatter(x=quad_v, y=bulkmod_y, mode='markers', marker=m3,
                       name='stencil'))


def eos1(mat: str, calc: str) -> None:
    pth = '/Users/ksb/scp_tmp/vauto/bulks/%s/%s/eos' % (calc, mat)
    print(os.listdir(pth))
    pe = [(x[0], x[2]) for x in zip(*parse_job(pth)) if x[-1] == '']
    assert pe, parse_job(pth)[-1]
    p, e = zip(*pe)

    v = [aseio.read(x + '/POSCAR').get_volume() for x in p]
    go.Figure([go.Scatter(x=v, y=e, mode='markers')]).show()


def vizone(mat: str, calc: str) -> None:
    q = '''SELECT name, volumes,energies,
                  expt_bm,bm,bulkmod_lstsq,expt_vol,lstsq_abc,cellvol,eng
           FROM bulks WHERE name=%s and calcname=%s'''
    name, v, e, expt_bm, bulkmod, bm_ls, expt_vol, \
        abc, vol, eng = dbgen.sqlselect(default.connect(), q, [mat, calc])[0]
    title, p2, p3, p4 = viz_plts(name, calc, v, e, bulkmod, vol,
                                 eng, abc, expt_vol, expt_bm, bm_ls)
    layout = go.Layout(hovermode='closest', showlegend=True, title=title,
                       xaxis=dict(title='Volume, À3'),
                       yaxis=dict(title='Energy, eV',))
    fig = go.Figure(data=[p2, p3, p4], layout=layout)
    plotly.offline.plot(fig, filename='temp0.html')


def viz(calc: str = 'ms2') -> None:
    q = '''SELECT name,volumes,energies,
                  expt_bm,bm,bulkmod_lstsq,expt_vol,lstsq_abc,cellvol,eng
           FROM bulks WHERE calcname=%s AND bm IS NOT NULL
           ORDER BY name'''
    res = dbgen.sqlselect(default.connect(), q, [calc])
    data, steps = [], []
    trues = [True, True, True]
    for i, (name, v, e, expt_bm, bulkmod, bm_ls, expt_vol,
            abc, vol, eng) in enumerate(res):

        title, plt2, plt3, plt4 = viz_plts(
            name, calc, v, e, bulkmod, vol, eng, abc, expt_vol,
            expt_bm, bm_ls)
        data.extend([plt2, plt3, plt4])

        args = [dict(visible=[False] * 3 * len(res)), {'title.text': title}]
        args[0]['visible'][3 * i:3 * i + 3] = trues  # type: ignore
        steps.append(dict(args=args, method='update'))

    layout = go.Layout(hovermode='closest',
                       xaxis=dict(title='Volume, À3'),
                       yaxis=dict(title='Energy, eV',),
                       showlegend=True, sliders=[dict(steps=steps)])
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename='temp0.html')

#############################################################################
###############
# FIT RELATED #
###############


def reset() -> None:
    q = '''TRUNCATE fitparams, fit'''
    dbgen.sqlexecute(default.connect(), q)
    os.system('rm -r /Users/ksb/scp_tmp/fitplt/*')


def mk_dataset(fx: str) -> fdata.Data:
    q = '''SELECT name, expt_ce,expt_bm,expt_vol,ab_ce,ab_bm,ab_vol, ce, volrat
           FROM bulks where calcname='{}' '''.format(fx)
    assert fx in ['ms2', '7500']
    datums = set()
    for (n, ce, bm, vol, abce_, abbm_, abvol_, ce_calc, volrat
         ) in dbgen.sqlselect(default.connect(), q):
        abce, abbm, abvol = map(lambda x: json.loads(   # type: ignore
            x) if x else '', [abce_, abbm_, abvol_])
        if abce and ce and ce_calc:
            datums.add(fdata.Datum(n, 'ce', abce[0], abce[1], float(ce), None))
        if abbm and bm:
            datums.add(fdata.Datum(n, 'bm', abbm[0], abbm[1], float(bm), None))
        if abvol and vol:
            datums.add(fdata.Datum(n, 'lc', abvol[0], abvol[1], float(vol),
                                   float(volrat)))

    return fdata.Data(datums)


def compare(fx: str, missing: bool = False) -> None:
    """
    Compare a functional to PBE/PBESOL/MS2/SCAN individually, only using the
    data points that they both have in common.
    """
    dsf = Dict[str, float]

    def get_data(f: str) -> T[dsf, dsf, dsf, List[int]]:
        ces, bms, lats = dict(), dict(), dict()
        nce, nbm, nlat = 0, 0, 0
        q = '''SELECT name, err_ce, expt_ce, relerr_bm, expt_bm, relerr_lat, expt_lat
               FROM bulks  WHERE calcname='{}' '''.format(f)
        for (n, ece, xce, ebm, xbm, elat, xlat) in dbgen.sqlselect(
                default.connect(), q):
            if xce is not None:
                nce += 1
                if ece is not None:
                    ces[n] = abs(float(ece))
                elif missing:
                    print(f, ' missing CE ', n)
            if xbm is not None:
                nbm += 1
                if ebm is not None:
                    bms[n] = abs(float(ebm))
                elif missing:
                    print(f, ' missing BM ', n)
            if xlat is not None:
                nlat += 1
                if elat is not None:
                    lats[n] = abs(float(elat))
                elif missing:
                    print(f, ' missing LC ', n)
        return ces, bms, lats, [nce, nbm, nlat]

    b1, b2, b3, ns = get_data(fx)
    base = [b1, b2, b3]
    for ref_ in ['pbe', 'pbesol', 'ms2', 'scan']:
        print('\n')
        if fx != ref_:
            r1, r2, r3, _ = get_data(ref_)
            ref = [r1, r2, r3]
            for (met, b, r, n) in zip(["c Energy", "(rel)bm ", "(rel)lat"],
                                      base, ref, ns):
                mats = set(b.keys()) & set(r.keys())  # Find overlapping mats
                nmat = len(mats)
                berr = round(sum(map(b.get, mats)) / nmat, 3)
                rerr = round(sum(map(r.get, mats)) / nmat, 3)
                print(met, '::', fx, ':', berr, ', ', ref_, ':', rerr,
                      '(miss. ', n - nmat, ')')


def errs(prop: str, rel: bool = True) -> None:
    '''Need to run `ms2fit` generator first'''
    assert prop in ['ce', 'bm', 'vol', 'lat']

    conn = default.connect()
    # dftval = dict(ce='ce',bm='bulkmod',l='lattice')[prop]
    relstr = 'rel' if rel else ''
    q1 = '''SELECT CONCAT('fit: ', CASE WHEN fit_id=-6532743347055789102 THEN 'ms2'
                                        WHEN fit_id=-505305133531848896 THEN 'pbesol'
                                        ELSE fit_id::Text END), x from fit'''
    q2 = '''SELECT CONCAT('fx: ',calcname), CONCAT('[',string_agg(
                    CONCAT('("',name,'",' ,{0}err_{1},')'),','),']')
            FROM bulks
            WHERE {0}err_{1} IS NOT NULL
            GROUP BY calcname'''.format(relstr, prop)
    fitx_ = dbgen.sqlselect(conn, q1)
    calcdata_ = dbgen.sqlselect(conn, q2)
    calcdata = []
    for k, v in calcdata_:
        mats, errs = zip(*ast.literal_eval(v))
        calcdata.append((k, mats, errs))

    fitx = [(k, np.array(json.loads(v))) for k, v in fitx_]
    key = dict(lat='lc', vol='lc')
    dataset = getattr(mk_dataset('ms2'), key.get(prop, prop))
    mats = [d.mat for d in dataset]

    fxdata = [go.Bar(name=k, x=mats, y=vals) for k, mats, vals in calcdata]

    fitdata = [go.Bar(name=n, x=mats,
                      y=[d.err(x, vol=prop == 'vol', rel=True)
                         for d in dataset]) for n, x in fitx]
    data = fxdata + fitdata

    fig = go.Figure(data=data, layout=go.Layout(barmode='group', title=prop,
                                                xaxis=dict(title='Material'), yaxis=dict(title=relstr + ' Error',)))
    plotly.offline.plot(fig, filename='temp0.html')


def maes(use_pareto: str = '') -> None:
    q = '''SELECT name, mae_ce,relmae_bm,relmae_lat FROM calc'''
    con = default.connect()
    fxfns, fxces, fxbms, fxlats = map(list, zip(*dbgen.sqlselect(con, q)))
    fxbumps = [-5 for _ in range(len(fxfns))]
    fxdatas = [fxces, fxbms, fxlats]

    dataset = mk_dataset('7500')

    fitfns, fitbumps, fitces, fitbms, fitlats = [], [], [], [], []
    q = '''SELECT fit_id, x, bump FROM fit'''
    for fitid, x_, bump in dbgen.sqlselect(con, q):
        fitfns.append(str(fitid))  # type: ignore
        fitbumps.append(bump or -5)
        x = json.loads(x_)
        fitces.append(dataset.mae(x, 'ce', rel=False))  # type: ignore
        fitbms.append(dataset.mae(x, 'bm', rel=True))  # type: ignore
        fitlats.append(dataset.mae(x, 'lc', rel=True))  # type: ignore
    fitdatas = [fitces, fitbms, fitlats]
    labs = [('' if x == 'CE' else "Relative")
            + ' MAE %s error' % x for x in ['CE', 'BM', 'Lattice']]
    pairs = [(0, 1), (0, 2)]

    data = []
    for x, y in pairs:
        dx, dy = copy.deepcopy(fitdatas[x]), copy.deepcopy(fitdatas[y])
        if use_pareto:
            breakpoint()
            dx, dyfnbumps = map(list, pareto(
                dx, list(zip(dy, copy.deepcopy(fitfns), copy.deepcopy(fitbumps)))))
            dy, dfitfns, dfitbumps = map(list, zip(*dyfnbumps))
        else:
            dfitfns = fitfns
            dfitbumps = fitbumps
        data.append(go.Scatter(
            x=dx + fxdatas[x], y=dy + fxdatas[y], mode='markers', text=dfitfns + fxfns,
            marker=dict(color=dfitbumps+fxbumps, colorscale='Viridis', size=14,
                        colorbar=dict(thickness=20)),
            hoverinfo='text'))

    fig = plotly.tools.make_subplots(rows=1, cols=2)
    for i, (x, y) in enumerate(pairs):
        # ;breakpoint()
        fig.append_trace(data[i], 1, i + 1)
        fig['layout']['xaxis%d' % (i + 1)].update(title=labs[x], zeroline=True)
        fig['layout']['yaxis%d' % (i + 1)].update(title=labs[y], zeroline=True,
                                                  rangemode='nonnegative')
    plotly.offline.plot(fig, filename='temp0.html')


def mk_plt(retry: bool = False, fit_id: str = None) -> None:
    condition = (" AND cast(fit_id as varchar) LIKE '%s%%%%'" %
                 fit_id) if fit_id else ''
    q = '''SELECT fit_id, x, cv, fitdata
            FROM fit JOIN calc   ON  calc=calc_id
            WHERE x IS NOT NULL''' + condition
    res = dbgen.sqlselect(default.connect(), q)
    print('plotting fx2d')
    for name_, x_, _, _ in res:
        name = str(abs(name_))
        d = '/Users/ksb/scp_tmp/fitplt/' + name
        if os.path.exists(os.path.join(d)) and not retry:
            print('skipping ', d)
            continue
        os.makedirs(d, exist_ok=True)

        x = fx.FromMatrix(np.array(json.loads(x_)), name=name)
        plotly.offline.plot(fx.plots([x]),
                            filename=os.path.join(d, 'fx2d.html'),
                            auto_open=False)
        plotly.offline.plot(fx.plots([x], long=True),
                            filename=os.path.join(d, 'fx2dlong.html'),
                            auto_open=False)
    print('other plots')
    for name_, x_, cv_, data_ in res:
        # DESERIALIZE
        name = str(abs(name_))
        cvunc, cvcon = json.loads(cv_) if cv_ else None, None
        x = fx.FromMatrix(np.array(json.loads(x_)), name=name)
        data = fdata.Data.from_list(json.loads(data_))

        # CONSTANTS
        d = '/Users/ksb/scp_tmp/fitplt/' + name
        if os.path.exists(os.path.join(d, 'fx.html')) and not retry:
            continue
        os.makedirs(d, exist_ok=True)
        mets = ['ce', 'bm', 'lc']
        scale = dict(ce=1, bm=100, lc=0.01)

        # HELPER FUNCS

        def n(nn: str) -> str:
            return os.path.join(d, nn + '.html')

        # PLOT RESIDUALS
        datalist = []
        for m in mets:
            mats, ys = [], []
            for dat in getattr(data, m):
                mats.append(dat.mat)
                ys.append(dat.err(x.x))
            datalist.append(go.Bar(x=mats, y=ys, name=m))

        layout = go.Layout(title='Residuals ', hovermode='closest',
                           xaxis=dict(title='Material'),
                           yaxis=dict(title='Error', ticklen=5, gridwidth=2))
        fig = go.Figure(data=datalist, layout=layout)
        plotly.offline.plot(fig, filename=n('resid'), auto_open=False)

        # PLOT FUNCTIONAL
        for deriv in ['fx', 'd2', 'ad2']:
            fig = x.plot3d(deriv=deriv, smax=4, amax=4)
            plotly.offline.plot(fig, filename=n(deriv),
                                auto_open=False)

        # PLOT CV DATA
        if cvcon and cvunc:
            ncv = len(cvunc['ce'])
            cvdata = [go.Bar(x=mets, marker_color='black',
                             y=[(cvunc[k][i] - cvcon[k][i]) / scale[k]
                                for k in mets])
                      for i in range(ncv)]
            layout = dict(title='CV MSE - ' + str(scale), barmode='group',
                          showlegend=False)
            cv = go.Figure(data=cvdata, layout=layout)
            plotly.offline.plot(cv, filename=n('cv'), auto_open=False)


def surf() -> None:
    q = '''SELECT xc, string_agg(CONCAT(mat, ' ', err),',')
           FROM surf GROUP BY xc'''
    mats = ['Rh', 'Pd', 'Pt', 'Cu', 'Co', 'Ru']
    data = []
    for xc, materr_ in dbgen.sqlselect(default.connect(), q):
        materr = [x.split() for x in materr_.split(',')]
        x, y = zip(*[(mats.index(mat), err) for mat, err in materr])
        data.append(go.Scatter(x=x, y=y, name=xc))

    layout = go.Layout(hovermode='closest',
                       xaxis=dict(title='Material'),
                       yaxis=dict(title='Site preference, eV',),)
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename='tmp', auto_open=False)


####################################################################
######################
# STATE RELATED #
######################
def sub(opt: str = '') -> None:
    '''Looks at morejobs and prints bash cmds to submit more.'''
    s = 'python -m functionals.CLI.submit '
    if opt:
        q = '''SELECT name, calcname FROM bulks WHERE success
               AND (optsuccess IS NULL or NOT optsuccess)'''
        for mat, calc in dbgen.sqlselect(default.connect(), q):
            print(s + '--mat={} --xc={}'.format(mat, calc))
    else:
        q = '''SELECT name, calcname, strains FROM bulks WHERE morejobs=1'''
        res = dbgen.sqlselect(default.connect(), q)
        pstr = s + '--mat={} --xc={} --strain="{} {}"'
        for mat, calc, strains in sorted(res):
            mstrain = json.loads(strains)[-1]
            print(pstr.format(mat, calc, mstrain, mstrain + 10))


def lattice() -> None:
    '''Writes opt_vol.txt by looking at current DB state.'''
    xtra = ','.join(['expt_%s IS NOT NULL' % x for x in ['ce', 'bm', 'lat']])
    q = """SELECT name, calcname, cellvol, lat, {} FROM bulks
        WHERE success AND vol IS NOT NULL""".format(xtra)
    with open(volcsv, 'w') as f:
        w = csv.writer(f)
        for row in sorted(dbgen.sqlselect(default.connect(), q)):
            w.writerow(row)


def expt() -> None:
    '''Aggregates experimental data'''
    root = datapth + '%s.csv'
    kjmol_to_ev = units.kJ / units.mol
    # mRyd_to_eV = 0.013605691455111256
    # Parse data from individual csvs
    # -------------------------------
    # with open(root % 'atom_eform_coh', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     elem_form = {k: float(v) for k, v in r}

    with open(root % 'keld', 'r') as f:
        r = csv.reader(f)
        next(r)
        keld = {k: (mbfloat(e), mbfloat(b), float(l), mbfloat(m))
                for k, s, e, b, l, m in r}

    with open(root % 'scheff', 'r') as f:
        r = csv.reader(f)
        next(r)
        scheff = {row[0]: tuple(map(float, row[1:]))
                  for row in r}

    # with open(root % 'kubaschewski', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     kub = {k: float(v) for k, _, v in r}

    # with open(root % 'rungs_tran', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     tran = {k: (mbfloat(a), mbfloat(b), mbfloat(c)) for k, a, b, c in r}

    # with open(root % 'errorestimate_lejaeghere', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     lej = {k: (float(a), float(b), float(c)) for k, a, b, c in r}

    # with open(root % 'cohesive_guillermet', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     guil = {k: float(v) for k, v in r}

    # with open(root % 'glasser_coh', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     glas = {k: float(v) for k, _, v in r}

    # with open(root % 'hsesol_schimka', 'r') as f:
    #     r = csv.reader(f)
    #     next(r)
    #     schim = {k: (float(bm), mbfloat(ce)) for k, _, bm, ce in r}

    with open(root % 'sol58lp_54coh', 'r') as f:
        r = csv.reader(f)
        next(r)
        sol = {k: (mbfloat(m), mbfloat(l), mbfloat(d), mbfloat(c),
                   mbfloat(b), r, u)
               for k, _, _, m, l, _, d, c, b, r, u in r}

    dicts = [sol, keld, scheff  # kub, tran, lej, schim, glas, sol
             ]  # type: List[Dict[str, Any]]
    # Make list of mats
    Mat = collections.namedtuple(
        'Mat', ['ce', 'bm', 'lc', 'mag'])

    mats = {k: Mat([], [], [], []) for k in itertools.chain(
        *[x.keys() for x in dicts])}  # type: Dict[str, Mat]

    # for k, (ecoh, bm, vol) in lej.items():
    #     mats[k].ce.append(ecoh*kjmol_to_ev)
    #     mats[k].bm.append(bm)
    #     struct = matdata[k].struct
    #     if struct != 'hcp':
    #         ndict = dict(bcc=2, fcc=4, diamond=8, zincblende=8, rocksalt=8,
    #                      cesiumchloride=2)
    #         mats[k].lc.append(((ndict[struct]*vol)**(1/3)))
    #     else:
    #         mats[k].lc.append((vol*matdata[k].ca_rat*(3**0.5)/4)**(1/3))
    # for k, vs in tran.items():
    #     for key, val in zip(['lc', 'bm', 'ce'], vs):
    #         getattr(mats[k], key).append(val)

    # for k, ce in guil.items():
    #     mats[k].ce.append(mRyd_to_eV * ce)

    # for k, ce in glas.items():
    #     mats[k].ce.append(kjmol_to_ev * ce / 2)

    # for k, (bm, c) in schim.items():
    #     mats[k].ce.append(kjmol_to_ev * c if c else None)
    #     mats[k].bm.append(bm)
    for k, (c, b, l, m) in keld.items():
        mats[k].bm.append(b)
        mats[k].ce.append(c)
        mats[k].lc.append(l)
        mats[k].mag.append(m)

    for k, (m, lx, debye, c, b, corr, unit) in sol.items():
        mats[k].mag.append(m)
        mats[k].bm.append(b)
        mats[k].lc.append(lx)
        if c is not None:
            if unit == 'kcal/mol':
                solce = c * units.kcal / units.mol
            elif unit == 'kJ/mol':
                solce = c * kjmol_to_ev
            elif unit == 'eV':
                solce = c
            else:
                raise ValueError((k, c, unit))
            if corr == 'False':
                assert debye
                solce += (9. / 8.) * units.kB * debye
            mats[k].ce.append(solce)

    for k, (l1, l2, _, _, l3, b1, b2, _, _, b3, c1, c2, _, _, c3
            ) in scheff.items():
        for z in ['lc', 'ce', 'bm']:
            getattr(mats[k], z).clear()  # remove other datatpoints
        mats[k].lc.append(l3 + l1 - l2)
        mats[k].ce.append(-(c3 + c1 - c2))  # flip sign
        mats[k].bm.append(b3 + b1 - b2)

    # for k, ce in kub.items():
    #     if not mats[k].ce:  # ONLY IF WE HAVE NO OTHER OPTION
    #         elems = re.findall(r'[A - Z][^A-Z]*', k)
    #         assert all([e in elem_form for e in elems]), elems
    #         extra = sum(map(elem_form.get, elems))
    #         assert extra
    #         mats[k].ce.append(kjmol_to_ev * (ce+extra) / 2)

    missingguess = set(mats.keys()) - set(matdata.keys())
    assert not missingguess, sorted(missingguess)
    missingdata = set(matdata.keys()) - set(mats.keys())
    print("WARNING! Missing data for %d mats: %s" %
          (len(missingdata), missingdata))
    # Write mats to expt csv
    with open(datapth + 'expt.csv', 'w') as f:
        w = csv.writer(f)
        for k, mat in sorted(mats.items()):
            w.writerow([k] + list(map(safeAvg, mat)))  # type: ignore

#############################################################################


def eos(calc: str, mat: str, e2: str = '') -> None:
    '''Visualize the first EOS directly from suncat cluster'''
    fst = ["ssh", "ksb@suncatls1.slac.stanford.edu"]
    root = '/nfs/slac/g/suncatfs/ksb/beefjobs/bulks/%s/%s/eos%s/strain_*/' % (
        calc, mat, e2)

    def splt(z): return list(filter(None, z.decode().split('\n')))
    vols = splt(subp.check_output(fst + ['grep "volume of" %sOUTCAR' % root]))
    engs = splt(subp.check_output(fst + ['tail -n 1 %sOSZICAR' % root]))
    volumes = []  # type: List[float]
    for vstr in vols:
        v = float(vstr.split()[-1])
        if v not in volumes:
            volumes.append(v)
    energies = [float(v.split()[4]) for v in engs[1::2]]
    vols, engs = zip(*sorted(zip(volumes, energies)))
    eos = ase.eos.EquationOfState(volumes, energies)
    _, _, eosbm = eos.fit()
    # compare to 5 point stencil
    dx = vols[1] - vols[0]
    stencil = (-1 * engs[4] + 16 * engs[3] - 30 * engs[2] +
               16 * engs[1] - engs[0]) / (12 * dx**2)
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3

    bulkmod = stencil * vols[2] * ev_a3_to_gpa  # APPROXIMATE
    eosbm = eosbm / ase.units.kJ * 1.0e24
    err = abs(bulkmod - eosbm) / eosbm
    if err > 0.20:
        print("Irregular!")
    print('dx', dx, '\nnvols ', vols, '\nengs ', engs, '\nbms ', bulkmod,
          eosbm)

    eos.plot(show=True)

#############################################################################


def _main() -> None:
    if sys.argv[1] in globals():
        func = globals()[sys.argv[1]]
        args = [] if len(sys.argv) < 2 else sys.argv[2:]
        func(*args)
    else:
        print('%s not in %s' % (sys.argv[1], list(
                                filter(lambda x: str.isalpha(x[0]),
                                       globals().keys()))))


if __name__ == '__main__':
    _main()
