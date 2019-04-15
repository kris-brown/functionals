# External modules
from sys     import argv
from sys     import argv
from io      import BytesIO
from base64  import decodebytes
from os      import environ
from os.path import join
from json    import loads, dumps, load
from math    import log10
from numpy   import array # type: ignore
from plotly.offline import plot # type: ignore
from plotly.graph_objs import Scatter, Figure, Layout # type: ignore
# Internal Modules
from dbgen import ConnectInfo, sqlselect
from functionals import FromMatrix,Functional,PBE,BEEF,SCAN
from functionals.fit.fit import Fit

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
def viz(stordir : str, cxn : ConnectInfo = default) -> None:
    '''Look at the E vs V curve for a material by name'''
    def decode_base64(data:bytes)->bytes:
        """Decode base64, padding being optional."""
        missing_padding = len(data) % 4
        if missing_padding != 0:
            data += b'='* (4 - missing_padding)
        return decodebytes(data)

    def viz_(img:bytes)->None:
        from matplotlib import pyplot as plt # type: ignore
        import matplotlib.image as mpimg # type: ignore
        i2 = mpimg.imread(BytesIO(img), format = 'png')

        plt.imshow(i2, interpolation = 'nearest')
        plt.show()

    q = '''SELECT img
            FROM bulks
                JOIN job ON job = job_id
            WHERE stordir=%s AND img IS NOT NULL
            LIMIT 1'''

    conn  = cxn.connect()
    res   = sqlselect(conn,q,[stordir])
    if not res: print('no img found'); return None
    image = res[0][0]
    i1    = decode_base64(bytes(image,encoding='UTF-8')) # type: ignore
    viz_(i1)
################################################################################
###############
# FIT RELATED #
###############

def fxviz(name:str,stepnum : str = None)->None:
    '''Compare multiple enhancement factor plots to PBE and BEEF'''
    q = '''SELECT pth,opt,a13,msb,consts
            FROM fit
                JOIN calc C ON calc=C.calc_id
                JOIN beef B ON B.functional=C.functional
                JOIN fitparams ON fitparams = fitparams_id
           WHERE name=%s AND opt IS NOT NULL'''
    res = sqlselect(default.connect(),q,[name])
    if res:
        for pth,opt,a1,msb,cs in res:
            name = name + '__' + str(float(a1))
            step = int(stepnum) if stepnum else loads(opt)[2]
            with open(pth+'/result.json','r') as fi:
                x = array(load(fi)[2][step][0])
            fx = FromMatrix(x,float(a1),float(msb),name+'[%s]'%cs)
    Functional.plots([fx,PBE,BEEF,SCAN])
    import pdb;pdb.set_trace()
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

def cviols(*fitnames:str)->None:
    q = '''SELECT data FROM fit WHERE name = %s'''
    for fn in fitnames:
        d = sqlselect(default.connect(),q,[fn])
        import pdb;pdb.set_trace()

def fxtraj(*fitnames:str)->None:
    '''
    Visualize the cost function and constraint violations over the iterations
    of a single fit
    '''
    print('Only doing mse_ce')
    q = '''SELECT pth,opt FROM fit WHERE name = %s '''
    dx = 10
    annotations,data = [],[]
    for fitname in fitnames:
        pth,opt_ = sqlselect(default.connect(),q,[fitname])[0]
        fit = Fit.from_json(pth)
        with open(pth+'/result.json','r') as fi: steps = load(fi)[2]
        ox_,_,ox = steps[loads(opt_)[2]]
        oy = fit.midcosts([ox_]*5)[0]
        #costs = [fit.midcosts(dumps([step]*5)) for step in steps[2]]
        #xs,ys = zip(*[(c,(msec+0+0)/1) for i,(msec,_,_,rc,rb,rl,c) in enumerate(costs) if i%1==0])
        xs,ys = zip(*[(cviol,fit.midcosts([x]*5)[0])
                        for i,(x,_,cviol) in enumerate(steps)])
        data.extend([Scatter(x=xs,y=ys,name=fitname,mode='markers',
                             text = [str(x) for x in range(len(xs))],
                             hoverinfo = 'text')])
        annotations.extend([dict(x=xs[0], y=ys[0], xref='x',yref='y',text='beef'),
                            dict(x=ox,    y=oy   , xref ='x',yref='y',text='opt')])

    layout = Layout(title= 'Trajectory of constrained optimization', hovermode= 'closest',
                    xaxis= dict(title= 'Constraint Violation', ticklen= 5, zeroline= False, gridwidth= 2),
                    yaxis= dict(title= 'MSE CE error, eV', ticklen= 5, gridwidth= 2,),
                    annotations = annotations)

    fig = Figure(data=data,layout=layout)

    plot(fig,filename='temp0.html')
    import pdb;pdb.set_trace()

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

funcs = set(['viz','unconverged','fxfiz','fxtraj','r2tradeoff','remove','add'])

################################################################################
def main() -> None:
    if argv[1] in globals():
        func = globals()[argv[1]]
        args = [] if len(argv) < 2 else argv[2:]
        func(*args)
    else:
        print('%s not in %s'%(argv[1],funcs))
if __name__=='__main__':
    main()
