# External
from typing import (Any,
                    Set      as S,
                    Dict     as D,
                    List     as L,
                    Tuple    as T,
                    Optional as O,
                    Callable as C)

from numpy        import array,empty,vstack # type: ignore
from random       import shuffle

'''
Data Preprocessing for nonlinear fitting
'''

################################################################################

class Datum(object):
    '''Something to be fit: a CE/BM/LC'''
    def __init__(self, mat : str, kind : str, vec : L[float], offset : float, target : float)->None:
        assert kind in ['ce', 'bm', 'lc']
        self.mat    = mat
        self.kind   = kind
        self.vec    = vec
        self.offset = offset
        self.target = target

    def __eq__(self, other:object) -> bool:
        return False if not isinstance(other,Datum) else \
             vars(self) == vars(other)

    def __hash__(self) -> int:
        return hash((self.mat,self.kind))

    def err(self, x : array, vol : bool = False) -> float:
        '''Compute error for this data point: for lc, compute either vol or lattice constant'''
        y = array(self.vec) @ x + self.offset
        if vol or self.kind != 'lc': return y - self.target
        else:
            hcp     = 'hcp' in self.mat
            factor  = 1./1.4142 if hcp else 1
            new_y   = (max(y,0)*factor)**(1/3)
            new_tar = (self.target*factor)**(1/3)
            return new_y - new_tar

class Data(object):
    def __init__(self, data : S[Datum]) -> None:
        self.ce = {d.mat:d for d in data if d.kind == 'ce'}
        self.bm = {d.mat:d for d in data if d.kind == 'bm'}
        self.lc = {d.mat:d for d in data if d.kind == 'lc'}

    def __eq__(self, other : object) -> bool:
        return False if not isinstance(other,Data) else \
             all([getattr(self,x)==getattr(other,x) for x in ['ce','bm','lc']])

    def __len__(self)->int: return len(self.ce) + len(self.bm) + len(self.lc)

    @property
    def data(self)->S[Datum]:
        return set(self.ce.values()) | set(self.bm.values()) | set(self.lc.values())

    def mse(self,x:array,key:str)->float:
        '''Mean squared error'''
        assert key in ['ce','bm','lc','vol']
        if   key == 'ce': return sum(d.err(x)**2 for d in self.ce.values())/len(self.ce)
        elif key == 'bm': return sum(d.err(x)**2 for d in self.bm.values())/len(self.bm)
        elif key == 'lc': return sum(d.err(x)**2 for d in self.lc.values())/len(self.lc)
        elif key =='vol': return sum(d.err(x,vol=True)**2 for d in self.lc.values())/len(self.lc)
        else: raise ValueError()

    def split(self)->T['Data','Data']:
        '''Partition a dataset into two equal size datasets'''
        ldata   = list(self.data)
        n       = len(self)
        allinds = list(range(n))
        shuffle(allinds)
        traininds,testinds = allinds[:n//2],allinds[n//2:]
        return (Data(set([ldata[i] for i in traininds])),
                Data(set([ldata[i] for i in testinds])))

    def xy(self, ce : float, bm : float, lc : float)->T[array,array]:
        '''Matrices for fast computation of cost function, weighted by scale'''
        X,Y = empty((0,64)),[]
        scale = dict(ce=ce,bm=bm,lc=lc)
        for d in self.data:
            s = scale[d.kind]
            assert s >= 0
            if s > 0:
                X = vstack((X,array([d.vec])/s))
                Y.append((d.target - d.offset)/s)
        return array(X),array(Y)

    @classmethod
    def from_list(cls,xs:list)->'Data':
        return cls({Datum(**x) for x in xs})

    def to_list(self)->L[dict]:
        return [vars(d) for d in self.data]
