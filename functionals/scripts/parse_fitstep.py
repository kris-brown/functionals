from typing      import List,Tuple, Callable as C
from os          import environ,listdir
from os.path     import join
from json        import load
from re          import findall, compile
#############################################################################
def parse_fitstep()->Tuple[List[str],List[int],List[float],List[float]]:
    '''
    Parse FITPATH to add logs of fitting to the DB
    '''
    # Helpers
    #--------
    def concat_map(f:C,xs:list)->list:
        out = [] # type: list
        for x in map(f,xs): out+=x
        return out

    reg = compile(r'\|\s+(\d+)\s+' # first group (integer)
                  r'\|.+?\|.+?'    # skip two cols
                  r'\|\s+([+-]?\d+\.\d+e[+-]\d+)\s+' # obj func (scientific)
                  r'\|.+?\|.+?'    # skip two cols
                  r'\|\s+([+-]?\d+\.\d+e[+-]\d+)' # const violation  (scientific)
                  r'.+?[$\n]')  # ignore rest of line

    root = environ['FITPATH']

    # Function applied to each folder
    #--------------------------------
    def process(name : str) -> List[Tuple[str,int,float,float]]:
        pth = join(root,name)

        # Parse log CSV manually
        with open(join(pth,'log'), 'r') as f:
            groups = findall(reg,f.read())

        return [(name,n,o,c) for n,o,c in groups]

    output = concat_map(process,listdir(root))
    names,niters,objs,cviols = map(list,zip(*output))

    return names,niters,objs,cviols # type: ignore

if __name__=='__main__':
    parse_fitstep()
