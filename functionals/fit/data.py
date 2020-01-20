# External
from typing import (Any,
                    Set as S,
                    Dict as D,
                    List as L,
                    Tuple as T)
import random
import math
import numpy as np

'''Data Preprocessing for nonlinear fitting.'''

###############################################################################

mags = ['MnC_b1', 'MnO_b1', 'FeN_b1', 'Fe_bcc', 'CrC_b1', 'CrN_b1', 'Ni_fcc',
        'FeAl_b2', 'Co_hcp', 'MnS_b1', 'MnN_b1']


class Datum(object):
    '''Something to be fit: a CE/BM/LC'''

    def __init__(self, mat: str, kind: str, vec: L[float], offset: float,
                 target: float) -> None:
        assert kind in ['ce', 'bm', 'lc']
        self.mat = mat
        self.kind = kind
        self.vec = vec
        self.offset = offset
        self.target = target

    def __eq__(self, other: object) -> bool:
        return False if not isinstance(other, Datum) else \
            vars(self) == vars(other)

    def __str__(self) -> str:
        return '<%s.%s>' % (self.mat, self.kind)

    def __hash__(self) -> int:
        return hash((self.mat, self.kind))

    def err(self, x: np.ndarray, vol: bool = False) -> float:
        ''' Compute error for this data point:
           - for lc, compute either vol or lattice constant
        '''
        y = np.array(self.vec) @ x + self.offset
        if vol or self.kind != 'lc':
            return float(y - self.target)
        else:
            hcp = 'hcp' in self.mat
            factor = 1./1.4142 if hcp else 1
            new_y = (max(y, 0)*factor)**(1/3)
            new_tar = (self.target*factor)**(1/3)
            return float(new_y - new_tar)


class Data(object):
    def __init__(self, data: S[Datum], full: bool = True) -> None:
        self.ce = sorted([d for d in data if d.kind == 'ce'],
                         key=lambda x: x.mat)
        self.bm = sorted([d for d in data if d.kind == 'bm'],
                         key=lambda x: x.mat)
        self.lc = sorted([d for d in data if d.kind == 'lc'],
                         key=lambda x: x.mat)
        assert self.ce and self.bm and self.lc

    def __eq__(self, other: object) -> bool:
        return False if not isinstance(other, Data) else \
            all([getattr(self, x) == getattr(other, x)
                 for x in ['ce', 'bm', 'lc']])

    def __len__(self) -> int: return len(self.ce) + len(self.bm) + len(self.lc)

    @property
    def data(self) -> S[Datum]:
        return set(self.ce) | set(self.bm) | set(self.lc)

    def mae(self, x: np.ndarray, key: str) -> float:
        '''Mean average error'''
        assert key in ['ce', 'bm', 'lc', 'vol']
        if key in ['ce', 'bm', 'lc']:
            return sum(abs(d.err(x)) for d in
                       getattr(self, key))/len(getattr(self, key))
        elif key == 'vol':
            return sum(abs(d.err(x, vol=True)) for d in self.lc)/len(self.lc)
        else:
            raise ValueError()

    def rmse(self, x: np.ndarray, key: str) -> float:
        '''Root mean square error'''
        assert key in ['ce', 'bm', 'lc', 'vol']
        if key in ['ce', 'bm', 'lc']:
            return (sum(d.err(x)**2 for d in
                        getattr(self, key))/len(getattr(self, key)))**0.5
        elif key == 'vol':
            return (sum(abs(d.err(x, vol=True)) for d in self.lc)/len(self.lc)
                    )**(0.5)
        else:
            raise ValueError()

    def xy(self, ce: float, bm: float, lc: float) -> T[np.ndarray, np.ndarray]:
        '''Matrices for fast computation of cost func, weighted by scale.'''
        X, Y = np.empty((0, 64)), []
        scale = dict(ce=ce, bm=bm, lc=lc)
        for d in self.data:
            s = scale[d.kind]
            assert s >= 0
            if s > 0:
                X = np.vstack((X, np.array([d.vec])/s))
                Y.append((d.target - d.offset)/s)
        return X, np.array(Y)

    @classmethod
    def from_list(cls, xs: L[Any], full: bool = False) -> 'Data':
        return cls({Datum(**x) for x in xs}, full=full)

    def to_list(self) -> L[D[str, Any]]:
        return [vars(d) for d in self.data]

    def split(self, n: int) -> L[T['Data', 'Data']]:
        '''Return a list of Leave One Out splits of the data.'''
        random.seed(42)
        ce, bm, lc = [sorted(getattr(self, x), key=str) for x in
                      ['ce', 'bm', 'lc']]
        random.shuffle(ce)
        random.shuffle(bm)
        random.shuffle(lc)

        def chunks(lst: L[Datum]) -> L[L[Datum]]:
            m = math.floor(len(lst)/n)
            return [lst[i*m:(i+1)*m] for i in range(n)]

        ces, bms, lcs = map(chunks, [ce, bm, lc])
        datas = []

        for i in range(n):
            test = set(ces[i] + bms[i] + lcs[i])
            testd = Data(test, full=False)
            traind = Data(self.data - test, full=False)
            datas.append((traind, testd))
        return datas
