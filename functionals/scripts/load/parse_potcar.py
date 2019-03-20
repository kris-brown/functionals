from typing import List as L, Tuple as T
from re     import search

def parse_potcar(pth : str) -> T[L[int],L[str],L[bool],L[int],L[float],      \
                                 L[float],L[float],L[float],L[float],        \
                                 L[float],L[float],L[bool],L[bool],L[float], \
                                 L[float],L[float],L[float],L[float]]:

    '''Parse a POTCAR file for a VASP job, which has an entry for each atom'''

    def boolean(x:str)->bool: return 'T' in x

    keys = dict(titel   = str,
                lultra  = boolean,
                iunscr  = int,
                rpacor  = float,
                pomass  = float,
                zval    = float,
                rcore   = float,
                rwigs   = float,
                enmax   = float,
                enmin   = float,
                lcor    = boolean,
                lpaw    = boolean,
                eaug    = float,
                rmax    = float,
                raug    = float,
                rdep    = float,
                rdept   = float) # type: dict


    def parse(s:str)->dict:
        '''Parse a POTCAR for a single atom by looking for keywords in a particular order'''
        lines = filter(lambda x: bool(x),   # ignore empty lines
                       iter(s.split('\n')))

        todo  = list(keys.items())[1:] # we do TITEL manually...since it has spaces in it

        currline = next(lines)
        d        = {'titel' : currline.strip()}
        currline = next(lines)

        while todo:
            kwl,func = todo.pop(0)
            kw    = kwl.upper()
            found = False
            while not found:
                if kw in currline:
                    p = kw+r'\s+=\s+(.*?)[;|\s|$]'
                    p = kw+r'\s+=\s+([^\s]+?)(;|\s|\Z)'
                    reg = search(p,currline)
                    assert reg
                    d[kw] = func(reg.groups()[0])
                    found = True
                else:
                    currline = next(lines)
        return d

    with open(pth,'r') as f:
        potcars = [parse(s) for s in f.read().split('End of Dataset')[:-1]]

    return tuple([list(range(len(potcars))),                         # type: ignore
                  *map(list,zip(*[p.values() for p in potcars]))])

if __name__=='__main__':
    print(parse_potcar('/Users/ksb/scp_tmp/vauto/bulks/AgCl_b1/strain_-1/POTCAR'))
