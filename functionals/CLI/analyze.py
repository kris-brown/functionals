# External modules
from typing  import List as L, Dict as D
from sys     import argv
from io      import BytesIO
from base64  import decodebytes
from os      import environ, system
from os.path import join
from time import sleep
from tempfile import NamedTemporaryFile
from json    import loads, dumps, load
from numpy   import array, linspace # type: ignore
from plotly import tools # type: ignore
from plotly.offline import plot # type: ignore
from plotly.graph_objs import Scatter, Figure, Layout # type: ignore
# Internal Modules
from dbgen import ConnectInfo, sqlselect
from functionals import FromMatrix,Functional,PBE,BEEF,SCAN
from functionals.fit.fit import Traj

"""
A CLI interface for common queries to the DB

Functions:
    - viz
    - unconverged
    - fxviz
    - fxtraj
"""


################################################################################
db      = '/'+join(*__file__.split('/')[:-3], 'data/functionals.json')
default = ConnectInfo.from_file(db)
ev_a3_to_gpa  = (10**-9) * (1.602 * 10**-19) * (10**10)**3

###############
# DFT RELATED #
###############

def viz(calc : int = 1) -> None:
    q = '''SELECT name, volumes,energies,allvol,alleng,eosbm,
                  expt_bm,bulkmod,bulkmod_lstsq,expt_vol,lstsq_abc,volume,eng
           FROM bulks WHERE calc=%s AND bulkmod IS NOT NULL
                        -- AND name in ('Li_bcc','Na_bcc')
           ORDER BY name'''
    res = sqlselect(default.connect(),q,[int(calc)])
    data, steps = [], []
    for i,(name, v_, e_, av_, ae_, eosbm, expt_bm, bulkmod_, bm_ls, expt_vol,
           abc, vol_, eng_) in enumerate(res):

        av, ae, v, e = map(loads,[av_,ae_,v_ or '[]',e_ or '[]'])
        bulkmod, vol, eng = map(float,[bulkmod_, vol_, eng_])
        m, m2, m3 = [dict(color = 'rgb(%d, %d, %d)'%(r,g,b),size = s,) for r,g,b,s
                     in [(231,99,250,20),(99,231,250,5),(199,231,50,5)]]

        title = '''{}<br><b>EOS</b>: eng={:.2f} eV, vol={:.1f} Å³, bm = {:.1f} GPa<br><b>Expt</b>: vol = {:.1f} Å³, bm = {:.1f} GPa<br><b>Stencil</b> bm: {:.1f} GPa, <b>QuadFit</b> bm: {:.1f} GPa'''
        title = title.format(name,eng,vol,eosbm,expt_vol or 0,expt_bm or 0, bulkmod, bm_ls)

        a,b,c = loads(abc)
        quad_v = linspace(min(v), max(v))
        quad_y = [a*x**2 + b*x + c for x in quad_v]
        bulkmod_a = bulkmod / (2*vol*ev_a3_to_gpa)
        bulkmod_y = [bulkmod_a*(x-vol)**2+eng for x in quad_v]

        data.extend([Scatter(x = av, y = ae,mode='markers', name='all'),
                     Scatter(x = v,  y = e, mode='markers', marker=m, name='opt'),
                     Scatter(x = quad_v,  y= quad_y, mode='markers',marker=m2, name='quad'),
                     Scatter(x = quad_v,  y= bulkmod_y, mode='markers',marker=m3,name='stencil')])

        args=[dict(visible=[False] * 4*len(res)),{'title.text':title}]
        args[0]['visible'][4*i:4*i+4] = [True, True, True, True] # type: ignore
        steps.append(dict(args=args,method='update'))

    layout = Layout(hovermode= 'closest',
                    xaxis= dict(title= 'Volume'),
                    yaxis= dict(title= 'Energy, eV',),
                    showlegend=True,
                        sliders=[dict(steps=steps)],)
    fig = Figure(data=data,layout=layout)
    plot(fig,filename='temp0.html')

################################################################################
###############
# FIT RELATED #
###############

def mses(constr : str = 'true') -> None:
    q = '''SELECT name,mse_ce,mse_bm,mse_lat FROM fit JOIN fitparams ON fitparams_id=fitparams WHERE {} UNION
           SELECT data,mse_ce,mse_bm,mse_lat FROM calc JOIN functional ON functional=functional_id WHERE mse_ce IS NOT NULL
        '''.format(constr)
    fns,ces,bms,lats = map(list,zip(*sqlselect(default.connect(),q)))
    datas = [ces,bms,lats]
    labs  = ['MSE CE, eV','MSE BM error, GPa','MSE lattice error, A']
    pairs = [(0,1),(0,2),(1,2)]
    data = [Scatter(x = datas[x], y = datas[y],  mode = 'markers',text = fns, hoverinfo = 'text') for x,y in pairs]

    fig = tools.make_subplots(rows=1, cols=3)
    for i,(x,y) in enumerate(pairs):
        fig.append_trace(data[i],1,i+1)
        fig['layout']['xaxis%d'%(i+1)].update(title=labs[x],zeroline=True,range=[0,0.2])
        fig['layout']['yaxis%d'%(i+1)].update(title=labs[y],zeroline=True,rangemode='nonnegative')
    plot(fig,filename='temp0.html')


def dcosts() -> None:
    q = sqlselect(default.connect(),'SELECT name,decaycosts,a11,a12,a13,a14,a15 FROM fit JOIN calc ON calc=calc_id JOIN functional ON calc.functional=functional_id JOIN beef ON beef.functional=functional_id')
    data = [] # type: list
    for n,d_,a1,a2,a3,a4,a5 in q:
        d = loads(d_)
        data.append(Scatter(x=[a1,a2,a3,a4,a5],y=[dd-d[2] for dd in d],name=n))
    layout = Layout(title= 'Decay parameter optimization', hovermode= 'closest',
                    xaxis= dict(title= 'Decay parameter'),
                    yaxis= dict(title= 'Relative weighted average R2 of fit'))
    fig = Figure(data=data,layout=layout);     plot(fig,filename='temp0.html')

################################################################################
######################
# SUBMISSION RELATED #
######################

def add() -> None:
    '''Prints bash statements to create extra jobs at needed strains'''
    q = '''SELECT stordir,morejobs,strain_low,strain_hi
            FROM bulks JOIN calc ON calc=calc_id JOIN functional ON functional=functional_id
            WHERE morejobs IN (-1,1)'''

    sroot= '/nfs/slac/g/suncatfs/ksb/functionals/data/structures/%s.traj'

    cmd = 'python functionals/CLI/submit.py --src={0} --time=40 --xc={1} --lo={2} --hi={3}'
    for dir,sgn,lo,hi in sqlselect(default.connect(),q):
        if   sgn == -1: lo    = abs(int(lo)-5)
        elif sgn == 1:  hi,lo = int(hi)+5,abs(lo)
        folders  = dir.split('/')
        material = sroot%folders[-1]
        xc       = folders[-2]
        cmd_ = cmd.format(material,xc.upper(),lo,hi)
        print(cmd_)
################################################################################
def main() -> None:
    if argv[1] in globals():
        func = globals()[argv[1]]
        args = [] if len(argv) < 2 else argv[2:]
        func(*args)
    else:
        print('%s not in %s'%(argv[1],globals().keys()))
if __name__=='__main__':
    main()
