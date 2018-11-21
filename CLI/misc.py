from typing    import Dict as D, List as L
from sys       import argv
from csv       import DictWriter
from os        import environ
from random    import choice
from operator  import mul
from functools import reduce

'''
Scripts to do misc things
'''
Dict = D[str,L]
################################################################################
def product(xs:L[int])->int:
    return reduce(mul,xs,1)

spec = dict(
    basis        = [2,2,3,3,3,4,4,5,6,7,8],
    initfit      = [0,1],
    bound        = [0.5,1],
    maxiter      = [1000,1500,2000],
    constden     = [3,5,7],
    constconst   = ['.']*8 + ['$a', "lda|pos"],       # '.' = True, '$a' = False
    dataconst    = ['.']*8 + ['^[A-Z][a-z]?_', '$a'], # single elements
    nlconstconst = ['.']*8 + ['$a'],                  # '.' = True, '$a' = False
    bm_weight    = [0.1,1,10],
    lat_weight   = [0.1,1,10]
) # type: Dict


class Counter(object):
    '''Basically a global variable counter'''
    val = 0
    @classmethod
    def get(cls)->int:
        cls.val += 1
        return cls.val

def fits(n:int=100)->None:
    with open(environ['FITPATH'],'a') as f:
        w = DictWriter(f,fieldnames=['name']+list(spec.keys()))
        w.writeheader()
        randfits(w,n)
        varyfits(w)

def varyfits(w : DictWriter)->None:
    '''Vary each parameter individually'''
    defaults = {k:v[len(v)//2] for k,v in spec.items()}
    w.writerow({**defaults,**{'name':'default'}})
    for k,v in spec.items():
        for val in set(v): # remove duplicates
            if val != defaults[k]: # don't repeat the base case
                w.writerow({**defaults,**{'name' : Counter.get(), k : val}})

def randfits(w : DictWriter, n : int) -> None:
    '''Randomly select n possible inputs'''
    for _ in range(n):
        w.writerow({**{'name':Counter.get()},
                    **{k:choice(v) for k,v in spec.items()}})
    tot = product([len(set(v)) for v in spec.values()])
    print('Sampled %d/%d possible configs'%(n,tot))


#############################################################################
#############################################################################
#############################################################################
def main() -> None:
    func = globals()[argv[1]]
    args = [] if len(argv) < 2 else argv[2:]
    func(*args)

if __name__=='__main__':
    main()
