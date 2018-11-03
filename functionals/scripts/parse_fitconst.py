from typing      import List,Tuple, Callable as C
from os          import environ,listdir
from os.path     import join
from json        import load
#############################################################################
def parse_fitconst() ->Tuple[List[str],List[str]]:
    '''
    Parse FITPATH to add logs of fitting to the DB
    '''
    # Helpers
    #--------
    def concat_map(f:C,xs:list)->list:
        out = [] # type: list
        for x in map(f,xs): out+=x
        return out

    root = environ['FITPATH']

    # Function applied to each folder
    #--------------------------------
    def process(name : str) -> List[Tuple[str,str]]:
        pth = join(root,name)

        with open(join(pth,'params.json'),'r') as f:
            params = load(f)

        return [(name,c) for c in params['consts']]

    output       = concat_map(process,listdir(root))
    names,consts = map(list,zip(*output))

    return names,consts # type: ignore

if __name__=='__main__':
    parse_fitconst()
