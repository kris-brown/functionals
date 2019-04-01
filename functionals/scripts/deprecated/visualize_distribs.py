from os.path           import exists,join
from os                import environ 
from sys               import argv
from json              import load
from numpy             import array # type: ignore
from numpy.random      import choice # type: ignore
from plotly            import tools # type: ignore
from plotly.graph_objs import Histogram as Hist,Layout,Figure # type: ignore
from plotly.offline    import plot # type: ignore

'''
Randomly sample s or a weighted by density to make a histogram
'''
################################################################################
def main()->None:
    # Params:
    N = 1000

    nbins = [100,1000]
    maxv  = [5,  20]
    # Init
    #-----
    names = ['s','alpha']
    rpt = lambda x: [x() for _ in range(N)]
    # Load input
    #-----------
    if len(argv) > 1:
        pth = argv[1]
        assert exists(pth), 'Point to a JSON output from analyize_density.py'
    else:
        pth = join(environ['HOME'],'scp_tmp/out.json')

    with open(pth,'r') as fi:
        dists = list(map(array,load(fi)))

    # Randomly sample with weighting, filter all below a threshold
    #----------------------------------------------------------
    #dists = [rpt(lambda:choice(a=mat.flatten(),p=pden)) for mat in [s,a]] # type: ignore

    # Plotly minutiae
    #----------------
    fig = tools.make_subplots(rows=1, cols=2, specs=[[{}, {}]],
                              shared_xaxes=False, shared_yaxes=False)


    data   = [Hist(x=dists[i],name=names[i],
                   xbins=dict(start=0,end=maxv[i],size=maxv[i]/nbins[i]),
                   histnorm ='probability')
                for i in range(2)]

    fig.append_trace(data[0], 1, 1)
    fig.append_trace(data[1], 1, 2)

    layout = Layout(xaxis=dict(title='s'),xaxis2=dict(title='a'))
    plot(fig)

if __name__=='__main__':
    main()
