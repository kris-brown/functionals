from typing import Dict

from dbgen.core.parsing import parse_line

def parse_xc_gpaw(log:str) -> str:
    """
    Finds the exchange correlation functional, bearing in mind the default
    """

    parsed =  parse_line(log,'xc:',0)
    if parsed is None:
        return 'LDA' # default
    else:
        return parsed.split(':')[-1].strip()
