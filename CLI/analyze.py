# External modules
from sys     import argv
from sys     import argv
from io      import BytesIO
from base64  import decodestring
from json    import loads
from numpy   import array # type: ignore
from plotly.offline import plot # type: ignore
from plotly.graph_objs import Scatter, Figure # type: ignore
# Internal Modules
from dbgen import ConnectInfo, sqlselect, hash_
from functionals import FromMatrix,Functional,PBE,BEEF,SCAN

"""
A CLI interface for common queries to the DB

Functions:
    - viz
    - unconverged
    - fxviz
    - fxtraj
"""
################################################################################
default = ConnectInfo.from_file('/Users/ksb/Documents/JSON/functionals.json')

def viz(material : str, cxn : ConnectInfo=default) -> None:

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

def fxviz(*fitnames:str)->None:
    '''Compare multiple enhancement factor plots to PBE and BEEF'''
    fxs = [] # type: list
    for fn in fitnames:
        q = '''SELECT basis,result FROM fit WHERE name=%s AND result IS NOT NULL'''
        res = sqlselect(default.connect(),q,[fn])
        if res:
            n,result = res[0]
            fxs.append(FromMatrix(array(loads(result)).reshape((n,n)),fn))
    Functional.plots(fxs+[PBE,BEEF,SCAN])

def fxtraj(*fitnames:str)->None:
    '''Visualize the cost function and constraint violations for a single fit'''

    def str2int(lst:str)->int:
        return sum(array(list(map(ord,hash_(lst)))))**3+10000000000

    q = '''SELECT FS.niter,FS.cost,FS.c_viol
            FROM fit_step FS JOIN fit ON fit_step__fit=fit_id
            WHERE fit.name = %s'''
    data = [] # type: list
    for fitname in fitnames:
        outs  = sqlselect(default.connect(),q,[fitname])
        i,x,y = zip(*outs)
        it    = zip([x,y],[max(x) or 1,max(y) or 1])          # iterate through this
        x_,y_ = [[item/mx for item in lst] for lst, mx in it] # normalized
        color = '#'+hex(str2int(fitname))[-6:]
        data.extend([Scatter(x=i,y=x_,name='%s-Cost'%fitname,
                             line={'color':color,'dash':'solid'}),
                     Scatter(x=i,y=y_,name='%s-Constraints'%fitname,
                             line={'color':color,'dash':'dash'})])
    fig   = Figure(data=data)
    plot(fig,filename='temp0.html')

def fxbest(start:int,stop:int,order:str='beefdist',kind:str='viz')->None:
    assert kind in ['viz','traj']
    assert order in ['resid','beefdist']
    i,j = map(int,[start,stop])
    q = '''SELECT fit.name FROM fit WHERE resid IS NOT NULL ORDER BY %s'''
    sel = sqlselect(default.connect(),q,[order])
    names = [x[0] for x in sel][i:j]
    if kind == 'viz':    fxviz(*names)
    elif kind == 'traj': fxtraj(*names)



funcs = set(['viz','unconverged','fxfiz','fxtraj'])
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
