# External Modules
from typing import List,Tuple,Any,Callable as C, Optional as O
from abc import ABCMeta,abstractmethod
from numpy import array # type: ignore
import numpy as np      # type: ignore
from scipy.optimize import LinearConstraint # type: ignore
from scipy.optimize import NonlinearConstraint # type: ignore

# Internal Modules
from functionals.fit.utilities import flatten, LegendreProduct

# Type synonyms #
Binary = C[[float,float],float]
################################################################################

################
# Helper funcs #
################

def const(x:Any)->C[[float,float],Any]:
    """Take a value and make a binary constant function"""
    return lambda _,__:x

true = const(True)

################
# Main classes #
################

class Constraint(metaclass=ABCMeta):
    """
    kind = {'lt','eq','gt','ltgt'}
    """

    @abstractmethod
    def make(self,n:int)->List[LinearConstraint]:
        raise NotImplementedError

    def __init__(self,
                 name : str, desc : str,
                 kind : str = 'eq')->None:
        assert kind in ['lt','gt','eq']
        self.kind = kind; self.name = name; self.desc = desc

    def bounds(self,v:float)->Tuple[float,float]:
        """ >, <, and = are all expressed using upper and lower bounds """

        if self.kind == 'eq':
            return v, v
        elif self.kind == 'lt':
            return -np.inf, v
        elif self.kind == 'gt':
            return v, np.inf
        else:
            raise ValueError


class MergedConstraint(Constraint):
    def __init__(self,cs:List['Constraint'])->None:
        self.cs = cs
        name = '_'.join(c.name for c in cs)
        desc = '\n\n'.join(c.desc for c in cs)
        super().__init__(name,desc,'eq')

    def make(self,n:int)->List[LinearConstraint]:
        """Prepare constraints for non-linear solver"""
        n2 = n**2
        lb,A,ub = np.empty((0,)),np.empty((0,n2)),np.empty((0,));
        for cons in flatten([c.make(n) for c in self.cs]):
            lb = np.concatenate((lb,cons.lb))
            A = np.vstack((A,array(cons.A)))
            ub = np.concatenate((ub,cons.ub))
        return LinearConstraint(A,lb,ub,keep_feasible=False)

    def linprog_matrices(self,n:int)->Tuple[array,array,array,array]:
        """Prepare constraints for linear programming solver"""
        n2    = n**2
        s1,s2 = (0,),((0,n2))
        A_eq,b_eq,A_ub,b_ub = np.empty(s2),np.empty(s1),np.empty(s2),np.empty(s1)

        for cons in flatten([c.make(n) for c in self.cs]):
            if cons.lb == cons.ub:
                A_eq = np.vstack((A_eq,array(cons.A)))
                b_eq = np.concatenate((b_eq,cons.ub))
            else:
                if cons.lb != -np.inf:
                    A_ub = np.vstack((A_ub, -1 * array(cons.A)))
                    b_ub = np.concatenate((b_ub, -1 * array(cons.lb)))
                if cons.ub != np.inf:
                    A_ub = np.vstack((A_ub, array(cons.A)))
                    b_ub = np.concatenate((b_ub, cons.ub))
        return A_eq,b_eq,A_ub,b_ub

class PointConstraint(Constraint):
    '''Express a constraint of Fx at some (s,⍺)'''
    def __init__(self,
                 name : str, desc : str,
                 s : float, alpha : float,val : float,
                 kind : str = 'eq')->None:
        self.s=s; self.alpha=alpha; self.val = val
        super().__init__(name,desc,kind)

    def make(self,n:int)->LinearConstraint:
        coefs = [[LegendreProduct(self.s,self.alpha,i,j) for j in range(n)] for i in range(n)]
        lo,hi = self.bounds(self.val)
        return [LinearConstraint([flatten(coefs)],[lo],[hi],keep_feasible=False)]

class FxConstraint(Constraint):
    """
    Express a constraint of the Fx(s,⍺) value (binary function) over some domain
    """
    def __init__(self,
                 name : str, desc : str,
                 ss : List[float], aa : List[float],f : Binary,
                 kind : str = 'eq')->None:
        self.ss=ss;self.aa=aa;self.f=f
        super().__init__(name,desc,kind)

    def make(self,n:int)->LinearConstraint:
        cons = []
        for s in self.ss:
            for a in self.aa:
                cons.append(PointConstraint('','',s,a,self.f(s,a),self.kind).make(n)[0])
        return cons



# Particular  Constraints
#-----------------------
cLDA = PointConstraint('LDA','LDA limit',s = 0, alpha = 1, val = 1, kind='eq')

default = [0.,0.1,0.5,1.,2.,3.,100.]
cLiebOx = FxConstraint('LiebOx','''The Lieb Oxford bound''',
    ss=default,aa=default,f=const(1.804),kind='lt')


cSCAN11 = FxConstraint('Scan11','''
Tighter L.O. bound (1.174) when ⍺=0 (Equation 11 in SCAN paper)''',
    ss=default,aa=[0],f=const(1.174),kind='lt')


cPos = FxConstraint('Positive','''Make Fx never negative''',
    ss=default,aa=default,f=const(0),kind='gt')

#######################
# NOT YET IMPLEMENTED #
#######################

# (at large s, say s=20): Fx(2s,0)/Fx(s,0) = 0.707
# but, at large s: s'(2s)==s'(s)
# Need to express as derivative condition???
#cSCAN13 = NonLinConstraint('Scan13','Scaling at large s (⍺=0) proportional to 1/sqrt(s) ')

# cGD     = [Constraint('GD',1e10,0,1.174)]#,LinConstraint(1e10,1e10,1.174)]
# cMS2    = Constraint('MS2',1e10,2,1.5)
