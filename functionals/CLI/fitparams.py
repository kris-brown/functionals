from typing    import Dict as D, List as L
from sys       import argv
from csv       import DictWriter
from os        import environ
from os.path   import join
from random    import choice
from operator  import mul
from functools import reduce
from itertools import product as prod

'''
Scripts to do misc things
'''

################################################################################
def product(xs : L[int]) -> int: return reduce(mul,xs,1)

spec = dict(
    constden   = [10],
    consts     = ['lda hnorm pos liebox scan11','lda hnorm pos', 'lda pos',''],
    reg        = [0,1],
    bm_weight  = [0],
    lat_weight = [0]
) # type: D[str,L]

def fits() -> None:
    with open(join(environ['FUNCTIONALS_ROOT'],'data/fitparams.csv'),'w') as f:
        w = DictWriter(f,fieldnames=list(spec.keys()))
        w.writeheader()
        allfits(w)# randfits(w,n); varyfits(w)

def allfits(w : DictWriter)->None:
    '''Full parameter search'''
    for combo in prod(*[set(x) for x in spec.values()]):
        kwargs = dict(zip(spec.keys(),combo))
        w.writerow(kwargs)

# class Counter(object):
#     '''Basically a global variable counter'''
#     val = 0
#     @classmethod
#     def get(cls)->int: cls.val += 1;  return cls.val

# def varyfits(w : DictWriter)->None:
#     '''Vary each parameter individually'''
#     defaults = {k:v[len(v)//2] for k,v in spec.items()}
#     w.writerow({**defaults,**{'name':'default'}})
#     for k,v in spec.items():
#         for val in set(v): # remove duplicates
#             if val != defaults[k]: # don't repeat the base case
#                 w.writerow({**defaults,**{'name' : Counter.get(), k : val}})
#
# def randfits(w : DictWriter, n : int) -> None:
#     '''Randomly select n possible inputs'''
#     for _ in range(n):
#         w.writerow({**{'name':Counter.get()},
#                     **{k:choice(v) for k,v in spec.items()}})
#     tot = product([len(set(v)) for v in spec.values()])
#     print('Sampled {}/{} possible configs'.format(n,tot))


if __name__=='__main__':
    fits()
