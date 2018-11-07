from typing      import List,Tuple
from os          import environ,listdir,mkdir
from os.path     import join,getmtime,exists
from json        import load,dumps
from datetime    import datetime # type: ignore
#############################################################################
def parse_fits()->Tuple[List[int],List[int],List[float]]:
    '''
    Parse FITPTH to add logs of fitting to the DB
    '''
    # Helpers
    #--------

    root = environ['FITPATH']

    # Function applied to each folder
    #--------------------------------
    def process(name : str) -> tuple:
        pth = join(root,name)

        # Get input parameters
        paramkeys = ['n','normc','ifit','bound','maxit','constden']

        with open(join(pth,'params.json'),'r') as f:
            params = load(f)

        base,norm,ifit,bound,maxiter,cden = map(params.get,paramkeys)

        # Get timestamp from folder
        time = datetime.fromtimestamp(getmtime(pth) + 8*3600) # add PST shift

        # Get result
        reskeys = ['fun','execution_time','x']
        with open(join(pth,'result.json'),'r') as f:
            result = load(f)

        resid,runtime,x = map(result.get,reskeys)

        return name,resid,runtime,dumps(x),time,base,norm,ifit,bound,maxiter,cden

    # MAIN #
    ########
    if not exists(root):
        mkdir(root)
    else:
        output = map(process,listdir(root))
        names,resids,runtimes,xs,times,bases,norms,ifits,bounds,maxiters,cdens,\
             = map(list,zip(*output))

    return names,resids,runtimes,xs,times,bases,norms,ifits,bounds,maxiters,cdens # type: ignore

if __name__=='__main__':
    parse_fits()
