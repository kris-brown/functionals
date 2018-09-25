# External modules
from typing import List, Tuple, Dict, Callable as C
from ast         import literal_eval
from math        import sqrt
from csv         import DictReader
from collections import defaultdict
import numpy as np        # type: ignore
from numpy import array   # type: ignore
# internal
from functionals.functional import beefcoeff
from functionals.utilities import corner
######################################################
class Data(object):
    """
    A target which should be minimized in a fitting problem
    """
    def __init__(self,x:array,target:array)->None:
        (r,c),(nt,) = x.shape,target.shape
        assert r==nt, '%d rows in X, %d elements in target'%(r,nt)
        self.n = nt; self.x = x; self.target = target

    @staticmethod
    def merge(ds:List['Data'])->'Data':
        """vertically stack a bunch of Data into one X,y"""
        x = np.vstack([d.x for d in ds])
        t = np.concatenate([d.target for d in ds])
        return Data(x,t)

    def A(self,n:int)->array:
        """Get subset of columns in the A matrix (from Ax=b) if our basis set
           has fewer than 8x8 elements."""
        sn = int(sqrt(n))
        assert sqrt(n) == sn, 'Expected %d to be a square number '%n
        return self.x[:,corner(sn)]

    def obj(self)->C[[array],float]:
        """
        Objective function for optimization: sum of squares of residual
        """

        def f(x:array)->float:
            A = self.A(len(x))
            resid = A @ x - self.target
            cost = resid @ resid
            return cost
        return f

    def jac(self)->C[[array],array]:
        """Jacobian of objective function (derived below in Einstein notation)

        Capital letters = vectors, matrices in [brackets]

        Fi  = residual                                                      --- (DIM: m x 1)
        Aij = basis function coefficients (cols) for each data point (row)  --- (DIM: m x n)
        Rj  = vector of weights corresponding to basis functions            --- (DIM: n x 1)
        Yi  = vector of targets for dot product of rows in Aij and Rj       --- (DIM: m x 1)
        δij = Kroenecker delta

        let: Fj = [A]ji * Ri - Yj
        d(Fj)/d(Rk) = [A]ji * d(Ri)/d(Rk) = [A]ji * δik = [A]jk

        d(FjFj)/d(Rk) = 2 * d(Fj)/d(Rk) * Fj
                      = 2 * [A]jk * Fj
        """
        def f(x:array)->array:
            A = self.A(len(x))
            residual = A @ x - self.target
            j = 2 * A.T @ residual
            return j
        return f

    def hes(self)->C[[array],array]:
        """
        Hessian of objective function: invariant w/r/t input vector
        To see this, take derivative of jac result w/r/t some Rm:

        d(FjFj)/d(Rk)d(Rm) = 2 * [A]jk * d(Fj)/d(Rm) =  2 * [A]jk * [A]jm
        """
        def f(x:array)->array:
            A = self.A(len(x))
            h = 2 * A.T @ A
            return h
        return f

def CohesiveData(pth:str)->Data:
    """
    Constructor for a Data object by pointing to CSV file with the columns:
        'Element','N_unit','keldAE','AtomE','BulkE','xAtom','xBulk'

    Could easily be made into a subclass of its own, if it were to need its own
    methods:

    def __init__(self,pth:str)->None:
        x,t = self.parseCSV(pth)
        super().__init__(x,t)
    """

    def parseCSV(pth:str)->Tuple[array,array]:
        """
        Parsing Outputs
            N x 64 matrices:
            Xbulk -  exchange contributions of optimized bulk structures per element
            Xatom -  exchange contributions of optimized isolated atoms per element

            vectors of length N:
            Ece   - cohesive energies per element from Keld
            Ebulk - total energies from optimized bulk structures per element
            Eatom - total energies from optimized isolated atoms per element
            Exbulk - exchange contributions to energies from optimized bulk structures per element
            Exatom - exchange contributions to energies from optimized isolated atoms per element
        """
        # Initialize lists/arrays + constants
        simpleCols = ['Element','N_unit','keldAE','AtomE','BulkE']
        simple = defaultdict(list) # type: Dict[str,list]
        XAtom,XBulk = [np.empty((0,8**2))]*2

        # Parse
        with open(pth, 'r') as f:
            reader = DictReader(f, delimiter=',', quotechar='"')
            for row in reader:

                # add the new row's values to the correct list
                for c in simpleCols:
                    simple[c].append(row[c])

                # 'flatten' the mxm matrices into vectors, then 'append'
                xAtom,xBulk = [array(literal_eval(row[y])).reshape((1,8**2))
                                    for y in ['xAtom','xBulk']]

                XAtom = np.vstack((XAtom,xAtom))
                XBulk = np.vstack((XBulk,xBulk))

        # Postprocess the resulting lists
        element,n_unit,target,Eatom,Ebulk = simple.values()

        elems   = element
        natoms  = array(n_unit).astype(int)
        Ece     = array(target).astype(float)
        Eatom   = array(Eatom).astype(float)
        Ebulk   = array(Ebulk).astype(float)
        Xbulk   = XBulk
        Xatom   = XAtom

        # Get real target
        #----------------
        Xbulk_norm = Xbulk / natoms[:,None]
        Ebulk_norm = Ebulk / natoms
        Exatom     = Xatom @ beefcoeff
        Exbulk     = Xbulk_norm @ beefcoeff
        dE         = (Ebulk_norm - Exbulk) - (Eatom - Exatom) # Non exchange contribution of cohesive energy (calculated by DFT, assumed valid)
        target     = Ece - dE                                 # 'True' difference E_x (bulk - atom)
        return Xbulk_norm,target

    x,t = parseCSV(pth)
    return Data(x,t)
