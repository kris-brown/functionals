import sys
import numpy as np
from typing import Tuple as T, Dict as D, Any, Union
import csv
import os
import plotly
import plotly.graph_objs as go
import json
import PyPDF2
from functionals.fit.functional import FromMatrix, SCAN, MS2
plotly.io.orca.config.executable = '/Applications/orca.app/Contents/MacOS/orca'
Arr = Union[np.ndarray]

mCaml = FromMatrix.frompath('msurf')

datapth = '/' + os.path.join(*__file__.split('/')[:-3], 'data/')

bulks = ['SOL54_coh', 'BM32', 'SOL58_lp']
gases = ['DBH24', 'RE42', 'S66x8']
surfs = ['CHEMI26', 'DISP15']
dsets = [bulks, gases, surfs]
alldsets = bulks + surfs + gases
fxs = [['mCAML', 'MS2', 'PBE', 'SCAN'],  # bulks
       ['mCAML', 'MS2', 'PBE', 'SCAN'],  # gases
       ['mCAML', 'PBE', 'SCAN', 'RPBE', 'MS2']]  # surfs

units = dict(SOL54_coh='eV', BM32='100 GPa', SOL58_lp='0.1 Å', DBH24='eV',
             CHEMI26='eV', RE42='eV', S66x8='0.1 eV', DISP15='eV')

scale = dict(BM32=100, S66x8=0.1, SOL58_lp=.1)
colors = dict(mCAML='red', MS2='green', SCAN='blue', PBE='black')
# MAE raw data
errs: D[str, D[str, float]] = dict(
    PBE=dict(SOL54_coh=.2, BM32=11.55, SOL58_lp=.051, DBH24=.33, DISP15=.49,
             CHEMI26=.27, RE42=.298,
             S66x8=.065, CHEMI26_me=-.17, DISP15_me=.49),
    mCAML=dict(SOL54_coh=.27, BM32=5.43, SOL58_lp=.032, DBH24=.5, DISP15=.16,
               CHEMI26=.22, RE42=.28,
               S66x8=.023, CHEMI26_me=-.17, DISP15_me=.06),
    MS2=dict(SOL54_coh=.23, BM32=4.34, SOL58_lp=.029, DBH24=.49, DISP15=.29,
             CHEMI26=.23, RE42=.42,
             S66x8=.052, CHEMI26_me=-.15, DISP15_me=.28),
    SCAN=dict(SOL54_coh=.21, BM32=4.39, SOL58_lp=.031, DISP15=.24, CHEMI26=.4,
              RE42=.338,
              S66x8=.0326, DBH24=.475, CHEMI26_me=-.4, DISP15_me=.04),
    RPBE=dict(CHEMI26=.17, DISP15=.72, RE42=.25, DBH24=.27, CHEMI26_me=.09,
              DISP15_me=.72))


def get_ab(dataset: str) -> T[Arr, Arr, Arr]:
    '''
    Read CSV file and apply normalization
    (to make different dataset errors comparable)
    Returns matrices for Ax=b
    '''

    # For spaghet, the MAE with real units is not the right number to normalize
    # b/c the raw results are given as percentages, so normalize these two
    # special cases by the unitless MAE.
    mae_msurf_lc = .0334
    mae_msurf_bm = .109

    assert dataset in alldsets
    A, b, c = [], [], []
    with open(datapth + 'msurf_contribs/' + dataset + '.csv',
              encoding='utf-8-sig') as f:
        r = csv.reader(f)
        for row in r:
            A.append(list(map(float, row[:64])))
            b.append(float(row[67]) - float(row[66]))
            c.append(float(row[68]) - float(row[66]))

    if dataset == 'SOL58_lp':
        err = mae_msurf_lc
    elif dataset == 'BM32':
        err = mae_msurf_bm
    else:
        err = errs['msurf'][dataset]

    prefactor = 1 / (len(b) * err)
    return (np.array(A) * prefactor,
            np.array(b) * prefactor,
            np.array(c) * prefactor)


def big_ab() -> T[Arr, Arr, Arr]:
    '''Big Ax=B matrices for all datasets to use for spaghetti'''
    A, b, c = np.empty((0, 64)), np.empty((0,)), np.empty((0,))
    for d in alldsets:
        A_, b_, c_ = get_ab(d)
        A = np.vstack((A, A_))
        b = np.concatenate((b, b_))
        c = np.concatenate((c, c_))

    with open(datapth + 'msurf_contribs/bigmatrix.json', 'w') as f:
        json.dump([A.tolist(), b.tolist(), c.tolist()], f)

    return A, b, c


def bivariate() -> None:
    fs = fxs[2]

    for mae in [True, False]:
        suf = '' if mae else '_me'
        A = 'A' if mae else ''
        emc = errs['mCAML']
        r = (emc['DISP15' + suf]**2 + emc['CHEMI26' + suf]**2
             )**0.5
        dispe = [errs[fx]['DISP15' + suf] for fx in fs]
        chem = [errs[fx]['CHEMI26' + suf] for fx in fs]
        fig = go.Scatter(x=dispe, y=chem, mode='markers',
                         marker_size=16 if mae else 12)
        arc = np.linspace(0, np.pi * (0.5 if mae else 2))
        curve = go.Scatter(x=r * np.cos(arc), y=r * np.sin(arc), mode='lines',
                           opacity=0.5)
        ax = dict(range=[0 if mae else -.75, .75 if mae else .8],
                  tickmode='linear', titlefont=dict(size=24 if mae else 18),
                  tickfont=dict(size=20 if mae else 16),
                  tick0=0, dtick=0.2 if mae else 0.3, zeroline=True,
                  gridcolor='grey' if mae else None, gridwidth=0.5,
                  zerolinecolor='grey', zerolinewidth=1)
        fig_ = go.Figure([curve, fig], layout=dict(
            showlegend=False,
            paper_bgcolor='white', plot_bgcolor='white',
            xaxis=dict(
                title='Dispersion M%sE (eV)' % A,
                side='bottom', **ax),
            yaxis=dict(title='Chemisorption M%sE (eV)' % A,
                       side='left', **ax)))
        wh = 600 if mae else 300
        fig_['layout'].update(width=wh, height=wh, autosize=False)

        off = dict(mCAML=(-.014, .05), MS2=(.055, .01), PBE=(.05, 0),
                   RPBE=(-.06, 0), SCAN=(0, .03)) if mae else dict(
                       mCAML=(-.35, 0), MS2=(0, -.17), PBE=(0.22, 0),
                       RPBE=(0, 0.17), SCAN=(0.2, -.17))
        for f, x, y in zip(fs, dispe, chem):
            dx, dy = off[f]
            fig_.add_annotation(
                text=f, x=x + dx, y=y + dy, showarrow=False,
                font=dict(size=18 if mae else 14, color="black"))
        fig_.write_image("data/figs/bivariate_M%sE.pdf" % A)

    # Make the figure-in-figure
    with open('data/figs/bivariate_MAE.pdf', 'rb') as f1:
        with open('data/figs/bivariate_ME.pdf', 'rb') as f2:
            with open('data/figs/bivariate.pdf', 'wb') as f3:
                p1 = PyPDF2.PdfFileReader(f1).getPage(0)
                p1.mergeTranslatedPage(PyPDF2.PdfFileReader(f2).getPage(0),
                                       200, 200)
                pdfWriter = PyPDF2.PdfFileWriter()
                pdfWriter.addPage(p1)
                pdfWriter.write(f3)


def err() -> None:
    '''
    Three subplots for the different types of datasets.
    For each, have a barplot where a column is a dataset and there are
    a group of bars, one for each functional.
    '''
    fig = plotly.subplots.make_subplots(rows=2, cols=1, subplot_titles=(
        "<b>Bulk properties</b>", "<b>Gas properties</b>"))

    for i in range(2):
        dset = dsets[i]
        for f in fxs[i]:
            fig.append_trace(go.Bar(
                name=f, x=['{}<br>({})'.format(d, units[d]) for d in dset],
                y=[errs[f][d] / scale.get(d, 1) for d in dset],
                marker_color=colors[f],
                legendgroup=1, showlegend=i == 1), i + 1, 1)

    ax = dict(gridcolor='grey', gridwidth=0.5, titlefont=dict(size=20),
              zerolinecolor='grey', zerolinewidth=1)
    fig.update_yaxes(title_text="MAE", row=1, col=1, **ax)
    fig.update_yaxes(title_text="MAE", row=2, col=1, **ax)
    fig['layout'].update(
        # title_font_size=50,
        font=dict(size=24),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)')

    fig.write_image("data/figs/err.pdf")


fsubx = '<i>f</i><sub>x</sub>'
its = '<i>s</i>'


def fx() -> None:
    srange = np.linspace(0, 4, 500)
    arange = np.linspace(0, 4, 500)
    fig = plotly.subplots.make_subplots(rows=2, cols=1)

    d1: D[Any, Any] = dict()  # xaxis=dict(title='s'))
    d2: D[Any, Any] = dict()  # xaxis=dict(title='α'))
    colors = ['red', 'blue', 'green']
    for c, fx in zip(colors, [mCaml, SCAN, MS2]):
        fig.append_trace(go.Scatter(
            x=srange, y=[fx.apply(s, 0) for s in srange],
            line=dict(color=c, dash='dash'),
            showlegend=False, **d1), 1, 1)
        fig.append_trace(go.Scatter(
            x=srange, y=[fx.apply(s, 1) for s in srange],
            line_color=c, showlegend=False, **d1), 1, 1)
        fig.append_trace(go.Scatter(
            y=[fx.apply(0, a) for a in arange], x=arange,
            line=dict(color=c, dash='dash'),
            showlegend=False, **d2), 2, 1)
        fig.append_trace(go.Scatter(
            y=[fx.apply(1, a) for a in arange], x=arange,
            line_color=c, showlegend=False, **d2), 2, 1)

    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='black', dash='dash'), name='α = 0'), 1, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='black', dash='solid'), name='α = 1'), 1, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='red', dash='solid'), name='mCAML'), 1, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='green', dash='solid'), name='MS2'), 1, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='blue', dash='solid'), name='SCAN'), 1, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='white', dash='solid'), name='<br>' * 1), 1, 1)

    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='black', dash='dash'), name=its + ' = 0'), 2, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='black', dash='solid'), name=its + ' = 1'), 2, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='red', dash='solid'), name='mCAML'), 2, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='green', dash='solid'), name='MS2'), 2, 1)
    fig.append_trace(go.Scatter(
        y=[0], x=[0], showlegend=True, visible='legendonly',
        line=dict(color='blue', dash='solid'), name='SCAN'), 2, 1)
    ax = dict(gridcolor='grey', gridwidth=0.5,
              zerolinecolor='grey', zerolinewidth=1)

    fig.update_yaxes(title_text=fsubx, row=1, col=1, **ax)
    fig.update_yaxes(title_text=fsubx, row=2, col=1, **ax)
    fig.update_xaxes(title_text=its, row=1, col=1, **ax)
    fig.update_xaxes(title_text="α", row=2, col=1, **ax)
    fig['layout'].update(
        font=dict(size=24),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)')
    fig.write_image("data/figs/fx.pdf")


def spaghet() -> None:
    fig = plotly.subplots.make_subplots(rows=3, cols=1, shared_xaxes=True,
                                        vertical_spacing=0.02)

    x = json.load(open(datapth + 'ensemble.json'))

    s = np.array(x[0])
    alphas = ['1', '0', '∞']
    fs = [np.array(x[i]) for i in [1, 3, 5]]
    ens = [np.array(x[i]) for i in [2, 4, 6]]
    for i, (alpha, f, en) in enumerate(zip(alphas, fs, ens)):
        fig.update_yaxes(
            title_text=fsubx + "(α=%s)" % (alpha), tick0=1, dtick=0.2,
            gridcolor='grey', row=i + 1, col=1)
        fig.update_xaxes(gridcolor='grey', row=i + 1, col=1)

        for e in en:
            fig.append_trace(go.Scatter(
                y=e, x=s, opacity=0.3, line_color='orange',
                showlegend=False), i + 1, 1)

        fig.append_trace(go.Scatter(
            x=s, y=f, line_color='black',
            showlegend=False), i + 1, 1)

    fig.update_xaxes(title_text=its, row=3, col=1)
    fig['layout'].update(font=dict(size=18),
                         paper_bgcolor='rgba(0,0,0,0)',
                         plot_bgcolor='rgba(0,0,0,0)')
    fig.write_image("data/figs/spaghet.pdf")


#############################################################################


def _main() -> None:

    if sys.argv[1] == 'all':
        funcs = ['spaghet', 'fx', 'err', 'bivariate']
    else:
        if sys.argv[1] in globals():
            funcs = [sys.argv[1]]
        else:
            print('%s not in %s' % (sys.argv[1], list(
                                    filter(lambda x: str.isalpha(x[0]),
                                           globals().keys()))))
            sys.exit()

    for f in funcs:
        func = globals()[f]
        args = [] if len(sys.argv) < 2 else sys.argv[2:]
        func(*args)
    print('DONE')


if __name__ == '__main__':
    _main()
