from dbgen.utils.parsing import parse_line

def parse_pw_gpaw(log:str) -> int:
    """
    Parse a planewave cutoff from GPAW logfile, bearing in mind its default values
    """

    # Parsing
    #-------
    rawmode = parse_line(log,'mode',0)
    if rawmode is None:
        raise ValueError('could not find what mode GPAW was run in')
    elif 'lcao' in rawmode:
        return 0 # non-planewave basis
    else:
        parsed  = parse_line(log,'ecut',0)
        if parsed is None:
            return 340 # gpaw default for PW mode
        else:
            raw_num = parsed.split()[-1].replace(',','')
            return int(round(float(raw_num)))
