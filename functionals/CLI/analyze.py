# External modules
from sys     import argv
from sys     import argv
from io      import BytesIO
from base64  import decodestring
from os      import environ
from os.path import join
from json    import loads, dumps
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
db      = join(environ['HOME'],'Documents/JSON/functionals_test.json')
default = ConnectInfo.from_file(db)

def viz(material : str, cxn : ConnectInfo = default) -> None:
    '''Look at the E vs V curve for a material by name'''
    def decode_base64(data:bytes)->bytes:
        """Decode base64, padding being optional."""
        missing_padding = len(data) % 4
        if missing_padding != 0:
            data += b'='* (4 - missing_padding)
        return decodestring(data)

    def viz_(img:bytes)->None:
        from matplotlib import pyplot as plt # type: ignore
        import matplotlib.image as mpimg # type: ignore
        i2 = mpimg.imread(BytesIO(img), format = 'png')

        plt.imshow(i2, interpolation = 'nearest')
        plt.show()

    q = '''SELECT img
            FROM functionals.expt
                JOIN species ON expt__species = species_id
            WHERE species.nickname=%s
            LIMIT 1'''

    conn  = cxn.connect()
    image = sqlselect(conn,q,[material])[0][0]
    i1    = decode_base64(bytes(image,encoding='UTF-8')) # type: ignore
    viz_(i1)

def unconverged() -> None:
    q = '''SELECT species_comp__element,symbol
            FROM species_comp JOIN element ON species_comp__element = element_id
            WHERE species_comp__element NOT IN
            (SELECT reference__element FROM
                reference JOIN calc ON reference__calc = calc_id
            WHERE xc = 'mBEEF')
            GROUP BY species_comp__element'''

    conn  = default.connect()
    print('\n'.join(map(str,sqlselect(conn,q))))

def fxviz(*datas:str)->None:
    '''Compare multiple enhancement factor plots to PBE and BEEF'''
    fxs = [] # type: list
    for fn in datas:
        q = '''SELECT result,a1,msb FROM fit
               WHERE name=%s AND result IS NOT NULL AND decay=2'''
        res = sqlselect(default.connect(),q,[fn])
        if res:
            for result,a1,msb in res:
                name = fn + '__' + str(float(a1))
                fxs.append(FromMatrix(array(loads(result)),float(a1),float(msb),name))
    Functional.plots(fxs+[PBE,BEEF,SCAN])

def r2tradeoff()->None:
    q = '''SELECT name,r2_ce,r2_bm,r2_lat,c_viol FROM fit'''
    res = sqlselect(default.connect(),q)
    data =[Scatter(x=[sum(x[1:4])/3 for x in res],y=[x[4] for x in res],mode='markers',
                         text = [x[0] for x in res])]
    layout = Layout(title= 'R2 vs constraint violation tradeoff', hovermode= 'closest',
                    xaxis= dict(title= 'Avg R2 value', ticklen= 5, zeroline= False, gridwidth= 2),
                    yaxis= dict(title= 'Constraint violation', ticklen= 5, gridwidth= 2,),
                    )
    fig = Figure(data=data,layout=layout)

    plot(fig,filename='temp0.html')

def fxtraj(*fitnames:str)->None:
    '''
    Visualize the cost function and constraint violations over the iterations
    of a single fit
    '''

    q = '''SELECT steps,fitparams,calc,decay
           FROM fit WHERE name = %s AND decay = 2'''

    annotations,data = [],[]
    for fitname in fitnames:
        steps,fp,calc,decay = sqlselect(default.connect(),q,[fitname])[0]
        fit = Fit.from_db(db,fp,calc,decay)

        costs = [fit.costs(db,step,calc,decay) for step in map(dumps,loads(steps))]
        xs,ys = zip(*[(log10(c),(rc+rb+rl)/3) for rc,rb,rl,c in costs])

        data.extend([Scatter(x=xs[5:],y=ys[5:],name=fitname,mode='markers+text',
                             text = list(map(str, range(len(xs)))),
                             textposition = 'bottom center')])
        annotations.append(dict(x=xs[-1],y=ys[-1],xref='x',yref='y',text='x',))


    layout = Layout(title= 'Trajectory of constrained optimization', hovermode= 'closest',
                    xaxis= dict(title= 'Log10(Constraint Violation)', ticklen= 5, zeroline= False, gridwidth= 2),
                    yaxis= dict(title= 'Average R2 value', ticklen= 5, gridwidth= 2,),
                    annotations = annotations)

    fig = Figure(data=data,layout=layout)

    plot(fig,filename='temp0.html')



funcs = set(['viz','unconverged','fxfiz','fxtraj','r2tradeoff'])


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
