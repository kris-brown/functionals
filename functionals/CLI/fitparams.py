from typing    import Dict as D, List as L
from sys       import argv
from csv       import DictWriter
from os        import environ
from os.path   import join
from itertools import product as prod

'''
Scripts to do misc things
'''

################################################################################
root = __file__.split('functionals')[0]
spec = dict(
    consts    = ['lda hnorm pos liebox scan11','lda pos','hnorm pos','lda hnorm pos', 'lda pos',''],
    reg       = [0,0.01,0.1],
    ce_scale  = [0.01,0.1],
    bm_scale  = [0.5,5],
    lc_scale  = [0.1,0.5]
) # type: D[str,L]

def fits() -> None:
    with open(join(root,'functionals/data/fitparams.csv'),'w') as f:
        w = DictWriter(f,fieldnames=list(spec.keys()))
        w.writeheader()
        allfits(w)# randfits(w,n); varyfits(w)

def allfits(w : DictWriter)->None:
    '''Full parameter search'''
    for combo in prod(*[set(x) for x in spec.values()]):
        kwargs = dict(zip(spec.keys(),combo))
        w.writerow(kwargs)

if __name__=='__main__':
    fits()
