# External modules
from typing  import List as L, Dict as D
from sys     import argv
from io      import BytesIO
from base64  import decodebytes
from os      import environ
from os.path import join
from json    import loads, dumps, load
from numpy   import array # type: ignore
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
db      = join(environ['FUNCTIONALS_ROOT'],'data/functionals.json')
default = ConnectInfo.from_file(db)

###############
# DFT RELATED #
###############

def viz(mat : str, calc : int = 1) -> None:
    q = 'SELECT volumes,energies,allvol,alleng,eosbm,expt_bm,bulkmod,lattice,expt_l FROM bulks JOIN job ON job=job_id WHERE name=%s and calc=%s'
    v_,e_,av_,ae_,eosbm,expt_bm,bulkmod,lattice,expt_l = sqlselect(default.connect(),q,[mat,calc])[0]
    av,ae,v,e = map(loads,[av_,ae_,v_,e_])
    m = dict(color = 'rgb(231, 99, 250)',size = 20,)
    d = dict(eosbm=eosbm,expt_bm=expt_bm,bulkmod=bulkmod,lattice=lattice,expt_l=expt_l)

    data  = [Scatter(x = av, y = ae,mode='markers',),
             Scatter(x = v,  y = e, mode='markers', marker=m)]

    title = mat + '  ' +'|'.join('' if v is None else '%s=%.2f'%(k,v)  for k,v in d.items())
    layout = Layout(title= title, hovermode= 'closest',
                    xaxis= dict(title= 'Volume'),
                    yaxis= dict(title= 'Energy, eV',),
                    showlegend=False)
    fig = Figure(data=data,layout=layout)

    plot(fig,filename='temp0.html')

################################################################################
###############
# FIT RELATED #
###############

def fxviz(name:str,stepnum : str = None)->None:
    '''Compare multiple enhancement factor plots to PBE and BEEF'''
    q = '''SELECT pth,result,a13,msb,consts
            FROM fit
                JOIN calc C ON calc=C.calc_id
                JOIN beef B ON B.functional=C.functional
                JOIN fitparams ON fitparams = fitparams_id
           WHERE name=%s'''
    res = sqlselect(default.connect(),q,[name])
    if res:
        for pth,xs,a1,msb,cs in res:
            name = name + '__' + str(float(a1))
            x = array(loads(xs)[2])
            fx = FromMatrix(x,float(a1),float(msb),name+'[%s]'%cs)
    Functional.plots([fx,PBE,BEEF,SCAN])
    #import pdb;pdb.set_trace()
def r2tradeoff()->None:
    q = '''SELECT name,r2_ce,r2_bm,r2_lat,c_viol FROM fit
            WHERE r2_ce{0}AND r2_bm{0}AND r2_lat{0}AND c_viol{0}'''.format(' IS NOT NULL ')
    res = sqlselect(default.connect(),q)
    data =[Scatter(x=[sum(x[1:4])/3 for x in res],y=[x[4] for x in res],mode='markers',
                         text = [x[0] for x in res])]
    layout = Layout(title= 'R2 vs constraint violation tradeoff', hovermode= 'closest',
                    xaxis= dict(title= 'Avg R2 value', ticklen= 5, zeroline= False, gridwidth= 2),
                    yaxis= dict(title= 'Constraint violation',),
                    )
    fig = Figure(data=data,layout=layout)

    plot(fig,filename='temp0.html')

def mses()->None:
    q = '''SELECT name,mse_ce,mse_bm,mse_lat FROM fit JOIN fitparams ON fitparams_id=fitparams UNION
           SELECT data,mse_ce,mse_bm,mse_lat FROM calc JOIN functional ON functional=functional_id WHERE mse_ce IS NOT NULL
        '''
    fns,ces,bms,lats = map(list,zip(*sqlselect(default.connect(),q)))
    datas = [ces,bms,lats]
    labs  = ['MSE CE, eV','MSE BM error, GPa','MSE lattice error, A']
    pairs = [(0,1),(0,2),(1,2)]
    data = [Scatter(x = datas[x], y = datas[y],  mode = 'markers',text = fns, hoverinfo = 'text') for x,y in pairs]

    fig = tools.make_subplots(rows=1, cols=3)
    for i,(x,y) in enumerate(pairs):
        fig.append_trace(data[i],1,i+1)
        fig['layout']['xaxis%d'%(i+1)].update(title=labs[x],zeroline=True,rangemode='nonnegative')
        fig['layout']['yaxis%d'%(i+1)].update(title=labs[y],zeroline=True,rangemode='nonnegative')
    plot(fig,filename='temp0.html')

def fxtraj(fitname : str, metric : str='cost') -> None:
    '''
    Visualize the cost function and constraint violations over the iterations
    of a single fit
    '''
    q = 'SELECT traj2 FROM fit WHERE name = %s'
    tstr = sqlselect(default.connect(),q,[fitname])[0][0]
    t = Traj.from_dict(loads(tstr))
    t.plot(metric)
    import pdb;pdb.set_trace()

def resid(fitname : str, metric : str = None) -> None:
    '''
    Visualize the cost function and constraint violations over the iterations
    of a single fit
    '''
    q = """SELECT traj2,CONCAT(pth,'/data.json'),ce_scale,bm_scale,lc_scale
           FROM fit JOIN fitparams ON fitparams_id=fitparams WHERE name = %s"""
    tstr,pth,ce,bm,lc = sqlselect(default.connect(),q,[fitname])[0]
    t = Traj.from_dict(loads(tstr),pth,2)
    met = [metric] if metric else None
    t.resid(met=met,scale=dict(ce=float(ce),bm=float(bm),lc=float(lc)))
    import pdb;pdb.set_trace()

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
def remove() -> None:
    '''Prints bash statements to remove directories with incomplete jobs'''
    atomq = 'SELECT stordir FROM atoms JOIN job ON job=job_id '\
            'WHERE not int_occupation or abs(true_mag - mag) > 0.05 '

    bulkq = 'SELECT stordir,incomplete FROM bulks JOIN job ON job=job_id '\
            'WHERE incomplete IS NOT NULL'

    for dir, in sqlselect(default.connect(),atomq):
        reldir = dir[25:]
        print('rm -rf '+reldir)
        #import pdb;pdb.set_trace()
    for dir,inc in sqlselect(default.connect(),bulkq):
        reldir = dir[25:]
        for i in map(int,inc.split()):
            print('rm -rf %s/strain_%d'%(reldir,i))

def add() -> None:
    '''Prints bash statements to create extra jobs at needed strains'''
    q = '''SELECT stordir,morejobs,strain_low,strain_hi
            FROM bulks JOIN job on job=job_id JOIN calc ON calc=calc_id JOIN functional ON functional=functional_id
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
# Error of fitted functional
# def besterr() -> None:
#     # Overall 'best'
#     best = array([0.9000268895332132, -0.5831983511873525, -0.5230678042284935, -0.06839485407616289, 0.0034552611415905048, -7.600333159159439e-05, 1.0705442345931335e-06, -1.1852935771847756e-08, 0.2835771188795937, 0.5210667876859505, 0.10134171361743657, -0.0032952119678765456, 5.1793960584492784e-05, -5.84916359952672e-07, 6.368179030973801e-09, -7.350101311747931e-11, 0.0315494983994424, 0.056031249847222867, -0.003119590227673151, 7.757592448868865e-05, -1.2211219001521302e-06, 1.4733208209958779e-08, -1.5932526149418185e-10, 1.6786540131206432e-12, -0.03895463463100185, -0.0001161277166103484, 2.0583447003345115e-05, -5.083937551897713e-07, 7.103784632017284e-09, -7.125548236945828e-11, 6.531811229342523e-13, -6.587646547481476e-15, -0.0020284224932872763, 1.7876789383377837e-05, -2.1782209662886135e-07, 3.082478394439253e-09, -3.0951971226390617e-11, 1.9659262742446466e-13, -1.2632129852520645e-15, 1.704045321242574e-17, -2.3934877774900652e-05, 2.122917591538513e-07, -5.016560962782947e-10, -2.4420350175009875e-11, 4.97039520150122e-13, -6.140406357139804e-15, 6.882186182099801e-17, -7.718526768936581e-19, 1.2835228838266299e-07, -6.195181082565845e-10, -1.8140010301737714e-11, 5.793235482140103e-13, -9.082298268357043e-15, 1.0275879330162299e-16, -1.0361503364890394e-18, 1.050616534279694e-20, 5.942187332047955e-09, -5.792415305747405e-11, 6.615003631873989e-13, -8.361699773265911e-15, 9.251606845820785e-17, -8.515824049703094e-19, 7.629426984934829e-21, -7.797107469268934e-23])
#     # best ce fx
#     best = array([373.70243900084125, 826.5757145613753, -572.425313475107, -2408.319018902571, -1042.2467678934072, 1409.8055037536103, 1082.3605618538802, 10.133162651191292, -1044.2453581882603, -2377.2205403823027, 1109.4892227296816, 5739.439137532255, 2114.4328967330844, -4150.567800987766, -3051.3657063067503, -79.78954935253904, 1278.5843411985802, 2886.77701975, -1504.9580821967747, -7184.030489264693, -2421.622179228543, 5481.357039349394, 3950.985171251446, 116.81766577875192, -1278.3822327309997, -3013.057536694557, 588.5329320630267, 5329.516439181601, 864.7852702488185, -5972.709942612807, -4032.65464546906, -206.24047918220884, 904.7066193531479, 2162.578938892508, -155.08221309176514, -3064.7659771616286, 335.9859971415498, 4916.893392665498, 3116.3018294451113, 193.7200202849139, -454.78869524636286, -1094.0277126410224, 72.7975747654176, 1425.474530698281, -559.6758373640013, -2956.1580478355136, -1808.0058723820334, -133.68268953578314, 260.93572341648587, 697.8573679659771, 481.8814760919519, 182.85624149666842, 988.9489227464614, 1668.4578339407046, 893.3252143865367, 80.05697034723678, 69.71389555827975, 215.12490078729715, 424.5741226989093, 493.0478304241885, 197.17393581873603, -93.24713960185152, -54.39824207172252, 19.902367498274504])
#     q = '''SELECT alloy,expt_ce,expt_bm,expt_vol,a_ce,a_bm,a_l,b_ce,b_bm,b_l FROM bulks WHERE a_ce IS NOT NULL and expt_ce IS NOT NULL'''
#     alloys = ['Metal','Hydride',"III",'IV',"V","VI","VII"]
#     errs = {k:[] for k in alloys} # type: D[str,L[float]]
#     tots = {k:0 for k in alloys}
#     for alloy,expt_ce,expt_bm,expt_vol,a_ce_,a_bm_,a_l_,b_ce_,b_bm_,b_l_ in sqlselect(default.connect(),q):
#         a_ce,a_bm,a_l = [array(loads(x)[2]) if x else None for x in [a_ce_,a_bm_,a_l_]]
#         b_ce,b_bm,b_l = [float(loads(x)[2]) if x else None for x in [b_ce_,b_bm_,b_l_]]
#
#         err = a_ce @ best + b_ce - float(expt_ce)
#         errs[alloy].append(err**2)
#         tots[alloy]+=1
#     rmses = {k:(sum(v)/tots[k])**(0.5) for k,v in errs.items()}
#     print(rmses)
#     import pdb;pdb.set_trace()
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
