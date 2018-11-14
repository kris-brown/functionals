from typing import Dict as D, List as L, Iterator as I
from sys import argv
from csv import DictWriter
from os import environ
from random import choice
from operator import mul
from functools import reduce
'''
Scripts to do misc things
'''
################################################################################
def product(xs:I[int])->int:
    return reduce(mul,xs,1)

def fits(n:int=100)->None:
    spec = dict(
        basis        = [2,2,3,3,3,4,4,5,6,7,8],
        initfit      = [0,1],
        bound        = [0.5,1],
        maxiter      = [1000,1500,2000],
        constden     = [3,5,7],
        constconst   = ['1']*8 + ['0', "const.const_name='lda'"],
        dataconst    = ['1']*8 + ['0', 'species.n_elems=1'],
        nlconstconst = ['1']*8 + ['0'],
        bm_weight    = [0.1,1,10],
        lat_weight   = [0.1,1,10]
    ) # type: D[str,L]
    with open(environ['FITPATH'],'w') as f:
        w = DictWriter(f,fieldnames=['name']+list(spec.keys()))
        w.writeheader()
        for i in range(n):
            w.writerow({**{'name':i},**{k:choice(v) for k,v in spec.items()}})

    print('Sampled %d/%d possible configs'%(n,product(map(len,spec.values()))))

def main() -> None:
    func = globals()[argv[1]]
    args = [] if len(argv) < 2 else argv[2:]
    func(*args)

if __name__=='__main__':
    main()
