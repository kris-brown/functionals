# External modules
from typing import Tuple as T, Any, Optional, Dict, List
import ast
import collections
import csv
import itertools
import os
import json
import sys
import numpy as np
import plotly
import plotly.graph_objs as go
import ase.units as units
import ase.io as aseio
# Internal Modules
import dbgen
import functionals.fit.functional as fx
import functionals.fit.data as fdata
from functionals.CLI.submit import matdata
from functionals.scripts.io.parse_job import parse_job

"""
A CLI interface for common queries to the DB
"""


#############################################################################
datapth = '/'+os.path.join(*__file__.split('/')[:-3], 'data/')
volcsv = datapth+'opt_vol.csv'
db = datapth+'functionals.json'
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
    #     print(xs); import pdb; pdb.set_trace();
    return sum(notnull)/len(notnull) if notnull else None

###############
# DFT RELATED #
###############


def viz_plts(name: str, calc: str,  v_: str, e_: str,
             bulkmod_: str, vol_: str, eng_: str, abc: str, eosbm: float,
             expt_vol: float, expt_bm: float,
             bm_ls: float) -> T[str, Any, Any, Any]:
    v, e = map(json.loads, [v_ or '[]', e_ or '[]'])
    bulkmod, vol, eng = [float(x) if x else 0. for x in [bulkmod_, vol_, eng_]]
    m, m2, m3 = [dict(color='rgb(%d, %d, %d)' % (r, g, b), size=s)
                 for r, g, b, s
                 in [(231, 99, 250, 20), (99, 231, 250, 5), (199, 231, 50, 5)]]

    title = '''{} {}<br><b>EOS</b>: vol={:.1f} Å³, bm = {:.1f}
                GPa<br><b>Expt</b>: vol = {:.1f} Å³, bm = {:.1f}
                GPa<br><b>Stencil</b> bm: {:.1f} GPa,
                <b>QuadFit</b> bm: {:.1f} GPa'''
    title = title.format(name, calc, vol, eosbm or 0, expt_vol or 0,
                         expt_bm or 0, bulkmod, bm_ls or 0)
    if abc:
        a, b, c = json.loads(abc)
        quad_v = np.linspace(min(v), max(v))
        quad_y = [a*x**2 + b*x + c for x in quad_v]
        bulkmod_a = bulkmod / (2*vol*ev_a3_to_gpa)
        bulkmod_y = [bulkmod_a*(x-vol)**2+eng for x in quad_v]
    else:
        quad_v = quad_y = bulkmod_y = []
    return (title,
            go.Scatter(x=v,  y=e, mode='markers', marker=m, name='opt'),
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

    v = [aseio.read(x+'/POSCAR').get_volume() for x in p]
    go.Figure([go.Scatter(x=v,  y=e, mode='markers')]).show()


def vizone(mat: str, calc: str) -> None:
    q = '''SELECT name, volumes,energies,eosbm,
                  expt_bm,bm,bulkmod_lstsq,expt_vol,lstsq_abc,vol,eng
           FROM bulks WHERE name=%s and calcname=%s'''
    name, v, e, eosbm, expt_bm, bulkmod, bm_ls, expt_vol, \
        abc, vol, eng = dbgen.sqlselect(default.connect(), q, [mat, calc])[0]
    title, p2, p3, p4 = viz_plts(name, calc, v, e, bulkmod, vol,
                                 eng, abc, eosbm, expt_vol, expt_bm, bm_ls)
    layout = go.Layout(hovermode='closest', showlegend=True, title=title,
                       xaxis=dict(title='Volume, À3'),
                       yaxis=dict(title='Energy, eV',))
    fig = go.Figure(data=[p2, p3, p4], layout=layout)
    plotly.offline.plot(fig, filename='temp0.html')


def viz(calc: str = 'ms2') -> None:
    q = '''SELECT name,volumes,energies,eosbm,
                  expt_bm,bulkmod,bulkmod_lstsq,expt_vol,lstsq_abc,volume,eng
           FROM bulks WHERE calc=%s AND bulkmod IS NOT NULL
           ORDER BY name'''
    res = dbgen.sqlselect(default.connect(), q, [calc])
    data, steps = [], []
    trues = [True, True, True, True]
    for i, (name, v, e, eosbm, expt_bm, bulkmod, bm_ls, expt_vol,
            abc, vol, eng) in enumerate(res):

        title, plt2, plt3, plt4 = viz_plts(
            name, calc, v, e, bulkmod, vol, eng, abc, eosbm, expt_vol,
            expt_bm, bm_ls)
        data.extend([plt2, plt3, plt4])

        args = [dict(visible=[False] * 4*len(res)), {'title.text': title}]
        args[0]['visible'][4*i:4*i+4] = trues  # type: ignore
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


def errs(prop: str) -> None:
    assert prop in ['ce', 'bm', 'vol', 'lat']
    abprop = 'vol' if prop == 'lat' else prop
    conn = default.connect()
    # dftval = dict(ce='ce',bm='bulkmod',l='lattice')[prop]

    q1 = '''SELECT CONCAT('fit: ',fit_id::Text), x from fit'''
    q2 = '''SELECT CONCAT('fx: ',calcname), CONCAT('[',string_agg(
                    CONCAT('("',name,'",' ,{0},',',expt_{0},')'),','),']')
            FROM bulks
            WHERE {0}-expt_{0} IS NOT NULL
            GROUP BY calcname'''.format(prop)
    q3 = '''SELECT name, ab_{0} FROM bulks WHERE calcname='ms2'
            AND ab_{0} IS NOT NULL AND expt_{1} IS NOT NULL
            '''.format(abprop, prop)
    fitx_ = dbgen.sqlselect(conn, q1)
    calcdata_ = dbgen.sqlselect(conn, q2)
    vecs_ = dbgen.sqlselect(conn, q3)

    fitx = [(k, np.array(json.loads(v))) for k, v in fitx_]
    calcdata_all = {k: ast.literal_eval(v) for k, v in calcdata_}
    calcdata = {k: {e[0]: e[1] for e in v} for k, v in calcdata_all.items()}
    exptvals = {e[0]: e[2]
                for e in set.union(*map(set, calcdata_all.values()))}

    vecs = {k: json.loads(v) for k, v in vecs_}
    mats = sorted([n for n, _ in vecs_] or set([k for v in calcdata.values()
                                                for k in v.keys()]))

    fxdata = [go.Bar(name=k, x=mats, y=[safe_sub(v.get(m), exptvals.get(m))
                                        for m in mats])
              for k, v in calcdata.items()]
    # import pdb;pdb.set_trace()
    fitdata = [go.Bar(name=n, x=mats,
                      y=[safe_sub(vecs[m][0]@x+vecs[m][1],
                                  exptvals.get(m)) for m in mats
                         ]) for n, x in fitx]
    data = fxdata + fitdata

    fig = go.Figure(data=data, layout=go.Layout(barmode='group'))
    # Change the bar mode
    plotly.offline.plot(fig, filename='temp0.html')


def maes() -> None:
    q = '''SELECT name,mae_ce,relmae_bm,relmae_lat,isfx FROM errs'''
    con = default.connect()
    fns, ces, bms, lats, fxs = map(list, zip(*dbgen.sqlselect(con, q)))
    datas = [ces, bms, lats]
    color = ['red' if x else 'blue' for x in fxs]
    labs = ['Relative MAE %s error' % x for x in ['CE', 'BM', 'Lattice']]
    pairs = [(0, 1), (0, 2)]  # ,(1,2)]
    data = [go.Scatter(x=datas[x], y=datas[y],  mode='markers', text=fns,
                       marker=dict(color=color),  hoverinfo='text')
            for x, y in pairs]

    fig = plotly.tools.make_subplots(rows=1, cols=2)
    for i, (x, y) in enumerate(pairs):
        # import pdb;pdb.set_trace()
        fig.append_trace(data[i], 1, i+1)
        fig['layout']['xaxis%d' % (i+1)].update(title=labs[x], zeroline=True)
        fig['layout']['yaxis%d' % (i+1)].update(title=labs[y], zeroline=True,
                                                rangemode='nonnegative')
    plotly.offline.plot(fig, filename='temp0.html')


def mk_plt(retry: bool = False) -> None:
    q = '''SELECT fit_id, x, cv, fitdata
            FROM fit JOIN calc   ON  calc=calc_id
            WHERE x IS NOT NULL'''
    res = dbgen.sqlselect(default.connect(), q)
    print('plotting fx2d')
    for name_, x_, _, _ in res:
        name = str(abs(name_))
        d = '/Users/ksb/scp_tmp/fitplt/'+name
        if os.path.exists(os.path.join(d)) and not retry:
            print('skipping ', d)
            continue
        os.makedirs(d, exist_ok=True)

        x = fx.FromMatrix(np.array(json.loads(x_)), name=name)
        plotly.offline.plot(fx.plots([x]),
                            filename=os.path.join(d, 'fx2d.html'),
                            auto_open=False)
    print('other plots')
    for name_, x_, cv_, data_ in res:
        # DESERIALIZE
        name = str(abs(name_))
        cvunc, cvcon = json.loads(cv_)
        x = fx.FromMatrix(np.array(json.loads(x_)), name=name)
        data = fdata.Data.from_list(json.loads(data_))

        # CONSTANTS
        d = '/Users/ksb/scp_tmp/fitplt/'+name
        if os.path.exists(os.path.join(d, 'fx.html')) and not retry:
            continue
        os.makedirs(d, exist_ok=True)
        mets = ['ce', 'bm', 'lc']
        scale = dict(ce=1, bm=100, lc=0.01)

        # HELPER FUNCS

        def n(nn: str) -> str:
            return os.path.join(d, nn+'.html')

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
        ncv = len(cvunc['ce'])
        cvdata = [go.Bar(x=mets, marker_color='black',
                         y=[(cvunc[k][i]-cvcon[k][i])/scale[k]
                            for k in mets])
                  for i in range(ncv)]
        layout = dict(title='CV MSE - '+str(scale), barmode='group',
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
            print(s+'--mat={} --xc={}'.format(mat, calc))
    else:
        q = '''SELECT name, calcname, strains FROM bulks WHERE morejobs=1'''
        res = dbgen.sqlselect(default.connect(), q)
        pstr = s+'--mat={} --xc={} --strain="{} {}"'
        for mat, calc, strains in sorted(res):
            mstrain = json.loads(strains)[-1]
            print(pstr.format(mat, calc, mstrain, mstrain+10))


def lattice() -> None:
    '''Writes opt_vol.txt by looking at current DB state.'''
    xtra = ','.join(['expt_%s IS NOT NULL' % x for x in ['ce', 'bm', 'lat']])
    q = """SELECT name, calcname, vol, lat, {} FROM bulks
        WHERE success AND vol IS NOT NULL""".format(xtra)
    with open(volcsv, 'w') as f:
        w = csv.writer(f)
        for row in sorted(dbgen.sqlselect(default.connect(), q)):
            w.writerow(row)


def expt() -> None:
    '''Aggregates experimental data'''
    root = datapth+'%s.csv'
    kjmol_to_ev = units.kJ/units.mol
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
                solce = c * units.kcal/units.mol
            elif unit == 'kJ/mol':
                solce = c * kjmol_to_ev
            elif unit == 'eV':
                solce = c
            else:
                raise ValueError((k, c, unit))
            if corr == 'False':
                assert debye
                solce += (9./8.)*units.kB*debye
            mats[k].ce.append(solce)

    for k, (l1, l2, l3, b1, b2, b3, c1, c2, c3) in scheff.items():
        for z in ['lc', 'ce', 'bm']:
            getattr(mats[k], z).clear()  # remove other datatpoints
        mats[k].lc.append(l3+l1-l2)
        mats[k].ce.append(-(c3+c1-c2))  # flip sign
        mats[k].bm.append(b3+b1-b2)

    # for k, ce in kub.items():
    #     if not mats[k].ce:  # ONLY IF WE HAVE NO OTHER OPTION
    #         elems = re.findall(r'[A-Z][^A-Z]*', k)
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
    with open(datapth+'expt.csv', 'w') as f:
        w = csv.writer(f)
        for k, mat in sorted(mats.items()):
            w.writerow([k] + list(map(safeAvg, mat)))  # type: ignore

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
