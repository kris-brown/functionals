

def all_mats(root: str) -> str:
    '''Get all materials that need CE, BM, and Lat data.

    This script needs to be in sync with parse_csv.'''
    from csv import reader

    ce, bml = set(), set()

    def add(m: str, c: str, b: str, v: str) -> None:
        if c:
            ce.add(m)
        if b or v:
            bml.add(m)

    with open(root % 'rungs_tran', 'r') as f:
        r = reader(f)
        next(r)
        for mat_, lat_, bm_, ce_ in r:
            add(mat_, ce_, bm_, lat_)

    with open(root % 'errorestimate_lejaeghere', 'r') as f:
        r = reader(f)
        next(r)
        for mat_, ce_, bm_, _ in r:
            if mat_ not in ce:
                add(mat_, ce_, bm_, '')

    with open(root % 'cohesive_guillermet', 'r') as f:
        r = reader(f)
        next(r)
        for mat_, ce_ in r:
            if mat_ not in ce:
                add(mat_, ce_, '', '')

    with open(root % 'sol58lp_54coh', 'r') as f:
        r = reader(f)
        next(r)
        for mat_, _, _, mag_, lp_, _, debye, ce_, bm_, corr, unit in r:
            add(mat_, ce_, bm_, lp_)

    if False:  # add keld data
        with open(root % 'misc', 'r') as f:
            r = reader(f)
            next(r)
            for mat_, prop, val, unit, _, _, _ in r:
                if prop == 'lattice':
                    add(mat_, '', '', val)
                elif prop == 'cohesive':
                    add(mat_, val, '', '')

    cc, bb = [' '.join(sorted(x)) for x in [ce, bml]]
    return cc + "|" + bb
