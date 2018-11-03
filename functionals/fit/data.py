# External modules
from typing      import List, Tuple, Dict, Callable as C
from ast         import literal_eval
from math        import sqrt
from csv         import DictReader
from collections import defaultdict
import numpy as np        # type: ignore
from numpy import array   # type: ignore

# Internal
from functionals.fit.functional import beefcoeff
from functionals.fit.utilities  import corner
from dbgen.support.misc         import ConnectInfo as Conn
from dbgen.core.sql             import select_dict
################################################################################
class Data(object):
    """
    A target which should be minimized in a fitting problem
    """
    def __init__(self,x:array,target:array)->None:
        (r,c),(nt,) = x.shape,target.shape
        assert r==nt, '%d rows in X, %d elements in target'%(r,nt)
        self.n = nt; self.x = x; self.target = target

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
            return np.zeros(len(x))
            # A = self.A(len(x))
            # h = 2 * A.T @ A
            # return h
        return f

    def linreg(self,n_:int)->Tuple[array,float]:
        """Solve Ax=b without constraints"""
        # Solve
        #-----
        n = n_**2

        a = np.linalg.lstsq(self.A(n),self.target,rcond=None)[0]
        resid = self.A(n) @ a - self.target
        return a, resid@resid

    @staticmethod
    def merge(ds:List['Data'])->'Data':
        """vertically stack a bunch of Data into one X,y"""
        x = np.vstack([d.x for d in ds])
        t = np.concatenate([d.target for d in ds])
        return Data(x,t)



class CohesiveData(Data):
    query = """
SELECT
    E.expt_id   AS uid,
    SP.nickname AS name,
    J.job_id    AS best_job_id,
    CA.coefs,
    SP.composition,
    CONCAT('{',
            GROUP_CONCAT(CONCAT(R.reference__element, ':', RJ.contribs)
                ORDER BY R.reference__element),
            '}') AS atomic_contribs,
    CONCAT('{',
            GROUP_CONCAT(CONCAT(R.reference__element, ':', RJ.energy)
                ORDER BY R.reference__element),
            '}') AS atomic_energies,
    J.contribs   AS bulk_contribs,
    J.energy     AS bulk_energy,
    SDE.value    AS target

FROM         expt                    E
        JOIN bulk_job                BJ  ON BJ.bulk_job__expt = E.expt_id
        JOIN job                     J   ON BJ.bulk_job__job  = J.job_id
        JOIN struct                  S   ON S.struct_id       = J.job__struct
        JOIN cell                    C   ON C.cell_id         = S.struct__cell
        JOIN species                 SP  ON E.expt__species   = SP.species_id
        JOIN species_comp            SC  ON SC.species_comp__species = SP.species_id
        JOIN reference               R   ON R.reference__element = SC.species_comp__element
                                            AND R.reference__calc = E.expt__calc
        JOIN job                     RJ  ON RJ.job_id     = R.reference__job
        JOIN calc                    CA  ON CA.calc_id    = J.job__calc
        JOIN species_dataset_element SDE ON SP.species_id = SDE.species_dataset_element__species
        JOIN (SELECT EE.expt_id,
                     MIN(ABS(EE.volume_pa - CC.volume / EE.n_atoms)) AS gap
              FROM    expt      EE
                 JOIN bulk_job  BBJJ ON BBJJ.bulk_job__expt = EE.expt_id
                 JOIN job       JJ ON BBJJ.bulk_job__job = JJ.job_id
                 JOIN struct    SS ON SS.struct_id = JJ.job__struct
                 JOIN cell      CC ON CC.cell_id = SS.struct__cell
              GROUP BY EE.expt_id
            ) AS X ON X.expt_id = E.expt_id
WHERE
    CA.xc = 'mBEEF'
    AND SDE.property = 'cohesive energy'
    AND ABS(E.volume_pa - C.volume / E.n_atoms) = X.gap
    AND eform IS NOT NULL
GROUP BY E.expt_id
ORDER BY nickname
"""
    def __init__(self,cxn : Conn) -> None:
        x,t = self._query(cxn)
        super().__init__(x,t)

    @classmethod
    def _query(cls,cxn : Conn) -> Tuple[array,array]:
        query_results = select_dict(cxn.connect(),cls.query)
        processed     = map(cls._process, query_results)
        xs,targets    = zip(*processed)
        import pdb;pdb.set_trace()
        return np.vstack(xs),array(targets)

    @staticmethod
    def _process(d:dict)->tuple:
        comp  = literal_eval(d['composition']) # type: Dict[int,int]

        ac_items        = literal_eval(d['atomic_contribs']).items()
        atomic_contribs = {i:array(xs).reshape((64,)) for i,xs in ac_items}
        ae_items        = literal_eval(d['atomic_energies']).items()
        atomic_energies = {i:float(e) for i,e in ae_items}

        e_atom = sum([atomic_energies[elem]*num for elem,num in comp.items()])
        x_atom = sum([atomic_contribs[elem]*num for elem,num in comp.items()])
        e_bulk = float(d['bulk_energy'])
        x_bulk = array(literal_eval(d['bulk_contribs'])).reshape((64,))
        dx     = x_atom - x_bulk
        coefs  = array(literal_eval(d['coefs'])) # allows us to convert xs -> energy

        # Use the BEEF coefficients from this particular calculator
        ex_atom = coefs @ x_atom
        ex_bulk = coefs @ x_bulk

        # Actual target is JUST the exchange component of formation energy
        nonx_e_atom = e_atom - ex_atom
        nonx_e_bulk = e_bulk - ex_bulk
        expt_co_eng = float(d['target']) # expt E atom - E bulk
        target = expt_co_eng  - (nonx_e_atom - nonx_e_bulk) # should be x_atom - x_bulk

        return (dx,target)
