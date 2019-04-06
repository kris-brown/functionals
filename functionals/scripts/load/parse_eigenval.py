from typing import Tuple as T

def parse_eigenval(fi:str)->T[bool,float]:
    '''Parse Vasp EIGENVAL - check if all bands have integer occupation'''
    magmom = 0.
    lines  = fi.split('\n')
    splitlines = filter(None,[x.split() for x in reversed(lines)])
    ints  = True
    for l in splitlines:
        if   len(l) == 3: vals = [float(l[-1])]              # spin unpolarized
        elif len(l) == 5: vals = [float(l[-2]),float(l[-1])] # spin polarized
        else: raise ValueError(l)
        nonints = [round(a,0) != round(a,2) for a in vals]
        magmom += vals[0] if (len(vals)==1) else abs(vals[0] - vals[1])

        if any(nonints): ints = False # flip flag if at any point this is true

        if float(l[0]) == 1.:     return ints,magmom # final band (reverse order)
    raise ValueError('Should not reach this part of code...')
