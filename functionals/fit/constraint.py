from typing import List as L, Dict as D,Tuple as T, Callable as C
from abc import ABCMeta,abstractmethod
from numpy import array,linspace,inf,empty,vstack,concatenate as concat # type: ignore
from numpy.linalg import norm # type: ignore
# Internal
from functionals.fit.utilities            import LegendreProduct
from functionals.scripts.fit.h_norm_const import intkernel,quad

Arrs       = L[array]
#####################################
def flatten(lol:list)->list: return [item for sublist in lol for item in sublist]

def h_norm_vec(a1 : float, msb : float) -> array:
    return [quad(intkernel(i,j,a1,msb), 0., inf)[0] for j in range(8) for i in range(8)]

class Constraint(object,metaclass=ABCMeta):
    def __init__(self,name:str,weight:float,kind:str)->None:
        self.name   = name
        self.weight = weight
        self.kind   = kind

        assert self.weight >0
        assert kind in ['lt','eq','gt']

    @abstractmethod
    def ab(self, decay:float, msb:float) -> T[array,array]: raise NotImplementedError
    @property
    @abstractmethod
    def card(self)->int: raise NotImplementedError

    @property
    def norm_weight(self)->float: return self.weight / self.card

    @classmethod
    def const_xy(cls,msb:float,decays:list,consts_:L['str'] = None) -> T[L[list],list,L[list],list]:
        '''
        Convert constraint dictionaries into matrices
        '''
        consts = [all_cons[x] for x in consts_] if consts_ is not None else list(all_cons.values())
        c_A_eqs,c_A_lts = [],[]

        for decay in decays:
            c_A_eq,c_b_eq,c_A_lt,c_b_lt = empty((0,64)),empty(0),empty((0,64)),empty(0)
            for const in consts:
                    A,b = const.ab(decay,msb)
                    if const.kind == 'eq':
                        c_A_eq = vstack((c_A_eq,A * const.norm_weight))
                        c_b_eq = concat((c_b_eq,b * const.norm_weight))
                    else:
                        c_A_lt = vstack((c_A_lt,A * const.norm_weight))
                        c_b_lt = concat((c_b_lt,b * const.norm_weight))

            c_A_lts.append([x.tolist() for x in c_A_lt])
            c_A_eqs.append([x.tolist() for x in c_A_eq])
        return c_A_eqs, c_b_eq, c_A_lts, c_b_lt

    def viol(self, x : array, a : float, m : float) -> float:
        '''Degree to which a BEEF vector violates the constraint'''
        A,b = self.ab(a,m)
        err = A@x - b
        if self.kind != 'eq': err = err * (err > 0) # Ax < b
        return sum(abs(err)) / self.card

class VC(Constraint):
    '''Vector constraint: one row ax = b, where a is a function of a1 and msb'''
    def __init__(self, name : str, weight : float, kind : str, fA : C[[float,float],float], val : float) -> None:
        self.fA = fA
        self.val = val
        super().__init__(name=name,weight=weight,kind=kind)

    @property
    def card(self)->int: return 1

    def ab(self, decay:float, msb:float) -> T[array,array]:
        return array([self.fA(decay,msb)]),array([self.val])

class LC(Constraint):

    def __init__(self,name:str,val:float,kind:str,weight:float=1.,
                 s:float=None,alpha:float=None,)->None:
        self.s        = s
        self.a        = alpha
        self.val      = val

        super().__init__(name=name,weight=weight,kind=kind)

    @property
    def grid(self)->L[float]: return linspace(0,5,20).tolist()
    @property
    def srange(self)->L[float]: return [float(self.s)] if self.s is not None else self.grid
    @property
    def arange(self)->L[float]: return [float(self.a)] if self.a is not None else self.grid
    @property
    def card(self)->int: return len(self.srange) * len(self.arange)

    def ab(self, decay : float, msb:float) -> T[array,array]:
        '''
        Turns a linear constraint into an Ax = b (or Ax < b) matrix
        '''
        sign   = -1 if self.kind == 'gt' else 1
        vals   = [self.val for _ in range(self.card)]
        c_A    = [flatten([[LegendreProduct(s,a,decay,msb,i,j) for j in range(8)]
                                                               for i in range(8)])
                                                            for s in self.srange
                                                            for a in self.arange]

        return array(c_A) * sign, array(vals) * sign

lda    = LC(name='lda',   weight = 1, s = 0,    alpha = 1,    kind = 'eq', val = 1.)
liebox = LC(name='liebox',weight = 1, s = None, alpha = None, kind = 'lt', val = 1.804)
scan11 = LC(name='scan11',weight = 1, s = None, alpha = 0.0,  kind = 'lt', val = 1.174)
pos    = LC(name='pos',   weight = 1, s = None, alpha = None, kind = 'gt', val = 0.0)

hnorm = VC(name='hnorm',weight=1.,kind='eq',fA=h_norm_vec,val= -.3125)

all_cons = dict(lda    = lda,
                liebox = liebox,
                scan11 = scan11,
                pos    = pos,
                hnorm  = hnorm) # type: D[str,Constraint]
