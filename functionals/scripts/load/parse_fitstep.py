from typing      import List as L, Tuple as T
from re          import findall

#############################################################################
def parse_fitstep(log:str)->T[L[str],L[int],L[float],L[float]]:
    '''
    Parse FITPATH to add logs of fitting to the DB
    '''
    # Helpers
    #--------

    reg = (r'\|\s+(\d+)\s+' # first group (integer)
           r'\|.+?\|.+?'    # skip two cols
           r'\|\s+([+-]?\d+\.\d+e[+-]\d+)\s+' # obj func (scientific)
           r'\|.+?\|.+?'    # skip two cols
           r'\|\s+([+-]?\d+\.\d+e[+-]\d+)' # const violation  (scientific)
           r'.+?[$\n]')  # ignore rest of line

    groups = findall(reg,log)
    niters,objs,cviols = map(list,zip(*groups))
    return niters,objs,cviols # type: ignore
