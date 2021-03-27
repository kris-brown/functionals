# External
from typing import (Any,
                    Set as S,
                    Dict as D,
                    List as L,
                    Tuple as T,
                    Optional as O)
import random
import math
import numpy as np
import plotly
import plotly.graph_objs as go

'''Data Preprocessing for nonlinear fitting.'''

###############################################################################
allmats = set(["AgBr", "AgCl", "AgF", "AlSb", "BaO", "BaSe", "CaS", "CaSe", "CdO", "CdS", "CdSe", "CdTe", "CoC", "CoN", "CrC", "CrN", "CsF", "CsI", "FeC", "FeN", "Ag", "AlAs", "AlN", "AlP", "Al", "Au", "BAs", "BN", "BP", "Ba", "C", "CaO", "Ca", "CoAl", "Cu", "FeAl", "Fe", "GaAs", "GaN", "GaP", "GaSb", "Ge", "HfC", "HfN", "InAs", "InP", "InSb", "Ir", "IrC", "IrN", "K", "KBr", "LaC", "LaN", "LiCl", "LiF",
               "LiH", "LiI", "Li", "MgO", "MgS", "MnC", "MnN", "MnO", "MnS", "Mo", "MoC", "MoN", "NaCl", "NaF", "Na", "NbC", "NbN", "Nb", "NiAl", "Ni", "NiC", "NiN", "OsC", "OsN", "Pd", "PdC", "PdN", "Pt", "PtC", "PtN", "Rb", "RbI", "Rh", "RhC", "RhN", "RuC", "RuN", "ScC", "ScN", "SeAs", "SiC", "Si", "Sn", "Sr", "Ta", "TaC", "TaN", "TiC", "TiN", "VC", "VN", "V", "W", "WC", "WN", "ZnS", "ZnSe", "ZnTe", "ZrC", "ZrN"])


class Datum(object):
    '''Something to be fit: a CE/BM/LC'''

    def __init__(self, mat: str, kind: str, vec: L[float], offset: float,
                 target: float, volrat: O[float]) -> None:
        assert kind in ['ce', 'bm', 'lc']
        assert (kind == 'lc') == bool(volrat)
        assert len(vec) == 64
        self.mat = mat
        self.kind = kind
        self.vec = vec
        self.offset = offset
        self.target = target
        self.volrat = volrat

    def __eq__(self, other: object) -> bool:
        return False if not isinstance(other, Datum) else \
            vars(self) == vars(other)

    def __str__(self) -> str:
        return '<%s.%s>' % (self.mat, self.kind)

    def __hash__(self) -> int:
        return hash((self.mat, self.kind))

    def err(self, x: np.ndarray, vol: bool = False, rel: bool = True) -> float:
        ''' Compute error for this data point:
           - for lc, compute either vol or lattice constant
        '''
        if len(self.vec) != 64 or len(x) != 64:
            breakpoint()
        y = (np.array(self.vec) @ x) + self.offset
        if vol or self.kind != 'lc':
            div = 1 / self.target if rel else 1
            return float(y - self.target) * div
        else:
            assert self.volrat
            new_y = (max(y, 0) / self.volrat)**(1 / 3)
            new_tar = (self.target / self.volrat)**(1 / 3)
            div = 1 / new_tar if rel else 1
            return float(new_y - new_tar) * div


class Data(object):
    def __init__(self, data: S[Datum], full: bool = True) -> None:
        self.ce = sorted([d for d in data if d.kind == 'ce'],
                         key=lambda x: x.mat)
        self.bm = sorted([d for d in data if d.kind == 'bm'],
                         key=lambda x: x.mat)
        self.lc = sorted([d for d in data if d.kind == 'lc'],
                         key=lambda x: x.mat)
        if full:
            assert self.ce and self.bm and self.lc

    def __eq__(self, other: object) -> bool:
        return False if not isinstance(other, Data) else \
            all([getattr(self, x) == getattr(other, x)
                 for x in ['ce', 'bm', 'lc']])

    def __len__(self) -> int:
        return len(self.ce) + len(self.bm) + len(self.lc)

    def remove(self, key: str, mat: str) -> bool:
        n = len(getattr(self, key))
        setattr(self, key, [x for x in getattr(self, key) if x.mat != mat])
        n2 = len(getattr(self, key))
        return n != n2  # successful removal

    @property
    def data(self) -> S[Datum]:
        return set(self.ce) | set(self.bm) | set(self.lc)

    def mae2(self, x: np.ndarray, key: str, rel: bool = True) -> float:
        '''This should give the same result as mae'''
        a, b = self._xy(getattr(self, key), rel)
        return float(np.sum(np.abs(a@x - b))) / len(b)

    def mae(self, x: np.ndarray, key: str, rel: bool = True) -> float:
        '''Mean average error'''
        assert key in ['ce', 'bm', 'lc', 'vol']
        if key in ['ce', 'bm', 'lc']:
            ln = len(getattr(self, key))
            if ln == 0:
                return 0.
            return sum(abs(d.err(x, rel=rel)) for d in
                       getattr(self, key)) / ln
        elif key == 'vol':
            if not self.lc:
                return 0.
            return sum(abs(d.err(x, rel=rel, vol=True))
                       for d in self.lc) / len(self.lc)
        else:
            raise ValueError()

    def rmse(self, x: np.ndarray, key: str) -> float:
        '''Root mean square error'''
        assert key in ['ce', 'bm', 'lc', 'vol']
        if key in ['ce', 'bm', 'lc']:
            return (sum(d.err(x)**2 for d in
                        getattr(self, key)) / len(getattr(self, key)))**0.5
        elif key == 'vol':
            return (sum(abs(d.err(x, vol=True)) for d in self.lc
                        ) / len(self.lc))**(0.5)
        else:
            raise ValueError()

    def xy(self, rel: O[bool] = None) -> T[L[np.ndarray], L[np.ndarray]]:
        '''Matrices for fast computation of cost func, weighted by scale.'''
        if rel is None:
            rels = [False, True, True]
        else:
            rels = [rel, rel, rel]
        return map(list, zip(*[self._xy(getattr(self, x), r)  # type: ignore
                               for x, r in zip(['ce', 'bm', 'lc'], rels)]))

    def _xy(self, ds: S[Datum], rel: bool) -> T[np.ndarray, np.ndarray]:
        X, Y = np.empty((0, 64)), []
        for d in ds:
            div = 1 / d.target if rel else 1
            X = np.vstack((X, np.array([d.vec]) * div))
            Y.append((d.target - d.offset) * div)
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
            m = math.floor(len(lst) / n)
            return [lst[i * m:(i + 1) * m] for i in range(n)]

        ces, bms, lcs = map(chunks, [ce, bm, lc])
        datas = []

        for i in range(n):
            test = set(ces[i] + bms[i] + lcs[i])
            testd = Data(test, full=False)
            traind = Data(self.data - test, full=False)
            datas.append((traind, testd))
        return datas

    def resid(self, key: str, x: np.ndarray, rel: bool = True) -> None:
        assert key in ['ce', 'bm', 'lc', 'vol']
        points = getattr(self, 'lc' if key == 'vol' else key)
        xs = [p.mat for p in points]
        y = [p.err(x, rel=rel, vol=key == 'vol') for p in points]
        layout = dict(title="%sError for %s" %
                      ('Relative ' if rel else '', key))
        fig = go.Figure(data=[go.Bar(name=key, x=xs, y=y)], layout=layout)
        plotly.offline.plot(fig, filename='temp0.html')
