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
A1 = 4.9479
def h_norm_vec(a1 : float=A1, msb : float=1) -> array:
    return [quad(intkernel(i,j,a1,msb), 0., inf)[0] for j in range(8) for i in range(8)]

class Constraint(object):
    '''A specific instance of a Constraint type (with a definite decay parameter and msb value)'''

    def __init__(self,name:str,A:array,b:array,eq:bool)->None:
        self.name = name; self.A = A; self.b = b; self.eq = eq

    def err(self,x:array)->float:
        '''
        The main purpose of a Constraint: to calculate the degree to which it
        is violated
        '''
        err = self.A@x - self.b
        if not self.eq: err = err * (err > 0) # Ax < b
        return sum(abs(err)) / len(self.b)

    def to_dict(self)->dict:
        return dict(name=self.name,A=self.A.tolist(),b=self.b.tolist(),eq=self.eq)

    @classmethod
    def from_dict(cls,d:dict)->'Constraint':
        return cls(name=d['name'],A=array(d['A']),b=array(d['b']),eq=d['eq'])

class ConstraintType(object,metaclass=ABCMeta):
    def __init__(self,name:str,kind:str)->None:
        self.name   = name
        self.kind   = kind
        assert kind in ['lt','eq','gt']

    def mkcon(self) -> Constraint:
        A, b = self.ab()
        return Constraint(self.name, A = A, b = b, eq = self.kind == 'eq')

    @abstractmethod
    def ab(self) -> T[array,array]: raise NotImplementedError

    @classmethod
    def const_xy(cls,consts_:L['str'] = None) -> T[L[list],list,L[list],list]:
        ''' Convert constraint dictionaries into matrices.'''
        consts = [all_cons[x] for x in consts_] if consts_ is not None else list(all_cons.values())

        c_A_eq,c_b_eq,c_A_lt,c_b_lt = empty((0,64)),empty(0),empty((0,64)),empty(0)
        for const in consts:
            A,b = const.ab()
            if const.kind == 'eq':
                c_A_eq = vstack((c_A_eq,A / len(b)))
                c_b_eq = concat((c_b_eq,b / len(b)))
            else:
                c_A_lt = vstack((c_A_lt,A / len(b)))
                c_b_lt = concat((c_b_lt,b / len(b)))

        return [x.tolist() for x in c_A_eq], c_b_eq, [x.tolist() for x in c_A_lt], c_b_lt

class VC(ConstraintType):
    '''Vector constraint: one row ax = b, where a is a function of a1 and msb'''
    def __init__(self, name : str, kind : str, vec : L[float], val : float) -> None:
        self.vec = vec
        self.val = val
        super().__init__(name=name,kind=kind)

    def ab(self,) -> T[array,array]:
        return array(self.vec),array([self.val])

class LC(ConstraintType):

    def __init__(self,name:str,val:float,kind:str,
                 s:float=None,alpha:float=None,)->None:
        self.s        = s
        self.a        = alpha
        self.val      = val

        super().__init__(name=name,kind=kind)

    @property
    def grid(self)->L[float]: return linspace(0,5,20).tolist() # CONSTRAINTS APPLIED FROM 0 TO 5 W/ 20 POINTS
    @property
    def srange(self)->L[float]: return [float(self.s)] if self.s is not None else self.grid
    @property
    def arange(self)->L[float]: return [float(self.a)] if self.a is not None else self.grid

    def ab(self) -> T[array,array]:
        '''
        Turns a linear constraint into an Ax = b (or Ax < b) matrix
        '''
        sign   = -1 if self.kind == 'gt' else 1
        c_A    = [flatten([[LegendreProduct(s=s,alpha=a,M=i,N=j,a1=A1,msb=1.) for j in range(8)]
                                                               for i in range(8)])
                                                            for s in self.srange
                                                            for a in self.arange]
        vals   = [self.val for _ in range(len(c_A))]

        return array(c_A) * sign, array(vals) * sign

lda    = LC(name='lda',    s = 0,    alpha = 1,    kind = 'eq', val = 1.)
liebox = LC(name='liebox', s = None, alpha = None, kind = 'lt', val = 1.804)
scan11 = LC(name='scan11', s = None, alpha = 0.0,  kind = 'lt', val = 1.174)
pos    = LC(name='pos',    s = None, alpha = None, kind = 'gt', val = 0.0)

hnorm = VC(name='hnorm',kind='eq',vec=h_norm_vec(),val= -.3125)

all_cons = dict(lda    = lda,
                liebox = liebox,
                scan11 = scan11,
                pos    = pos,
                hnorm  = hnorm) # type: D[str,ConstraintType]
