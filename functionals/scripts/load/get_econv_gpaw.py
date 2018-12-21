from typing import Optional as O
from dbgen.utils.parsing import parse_line
######################################################
def get_econv_gpaw(log : str) -> O[float]:
    """
    Electronic energy convergence
    """
    init = log.split('Initialize')[0]
    parsed = parse_line(init,'energy:',0)
    if parsed is None: # the first hit for 'energy:' is NOT in the input parameter section
        return 0.0005 # default
    else:
        raw    = parsed.split('energy:')[-1]
        pure   = raw.replace('}','')
        return float(pure)
