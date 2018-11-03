# External
from typing             import Tuple,Dict,List as L, Any
from ast                import literal_eval
from numpy              import vstack,array,logspace,inf # type: ignore
from json               import dump
from shutil             import copy
from os                 import mkdir
from os.path            import exists,join

# Internal
from dbgen.core.misc           import hash_
from dbgen.core.sql            import select_dict
from dbgen.support.misc        import ConnectInfo as Conn
from dbgen.core.lists          import flatten
from functionals.fit.utilities import LegendreProduct
###############################################################################
def safeMkdir(pth:str)->None:
    if not exists(pth): mkdir(pth)

class Fit(object):
    def __init__(self, cxn: Conn, size : int, bound: float = 0.1,
                 norm : float = 0.1, initfit : bool = True,
                 maxiter : int = 1000, gridden : int = 5 ) -> None:
        # Store inputs
        #-------------
        self.size    = size
        self.gridden = gridden

        # Query DB and store results
        #---------------------------
        q1   = 'SELECT * FROM cohesive_data WHERE %s'
        q2   = 'SELECT * FROM const WHERE %s'
        conn = cxn.connect()
        self.cohesive_data = select_dict(conn,q1,[1])
        self.consts        = select_dict(conn,q2,[1])


        self.args = dict(n     = size, bound = bound,  maxit    = maxiter,
                         normc = norm, ifit = initfit, constden = gridden,
                         consts = [d['name'] for d in self.consts])

        self.name = hash_(self.args)[:20]

    def data(self)->Tuple[L[L[float]],L[float]]:
        '''Returns A and b for an Ax = b fitting problem'''
        processed   = map(self._process, self.cohesive_data)
        xs,targets  = zip(*processed)
        return vstack(xs).tolist(), list(targets)

    def _process(self,d:dict)->tuple:
        comp  = literal_eval(d['composition']) # type: Dict[int,int]

        def mkarray(xs:list) -> array:
            return array(xs)[:self.size,:self.size].reshape((self.size**2))

        ac_items        = literal_eval(d['atomic_contribs']).items()
        atomic_contribs = {i:mkarray(xs) for i,xs in ac_items}
        ae_items        = literal_eval(d['atomic_energies']).items()
        atomic_energies = {i:float(e) for i,e in ae_items}

        e_atom = sum([atomic_energies[elem]*num for elem,num in comp.items()])
        x_atom = sum([atomic_contribs[elem]*num for elem,num in comp.items()])

        # allows us to convert xs -> energy
        coefs  = mkarray(array(literal_eval(d['coefs'])).reshape((8,8)))

        rat    = int(d['bulk_ratio'])
        e_bulk = float(d['bulk_energy']) / rat
        x_bulk = mkarray(literal_eval(d['bulk_contribs'])) / rat

        dx     = x_atom - x_bulk # THIS IS WHAT WE ARE FITTING

        # Use the BEEF coefficients from this particular calculator
        ex_atom = coefs @ x_atom
        ex_bulk = coefs @ x_bulk

        # Actual target is JUST the exchange component of formation energy
        nonx_e_atom = e_atom - ex_atom
        nonx_e_bulk = e_bulk - ex_bulk
        expt_co_eng = float(d['target']) # expt E atom - E bulk
        target = expt_co_eng  - (nonx_e_atom - nonx_e_bulk) # should be x_atom - x_bulk

        return (dx,target)

    def const(self)->Tuple[L[L[float]],L[float],L[float]]:
        '''Returns A, lb, and ub for linear constraints: lb <= Ax <= ub'''

        keys = ['s','alpha','val','kind']
        grid = logspace(-2,2,self.gridden).tolist()
        coefs,lo,hi = [],[],[] # type: Tuple[list,list,list]
        for s_,a_,val_,kind in [map(c.get,keys) for c in self.consts]:
            val    = float(val_)
            srange = [float(s_)] if s_ else grid
            arange = [float(a_)] if a_ else grid
            if   kind == 'eq': lo_,hi_ = val, val
            elif kind == 'lt': lo_,hi_ = -inf, val
            elif kind == 'gt': lo_,hi_ = val, inf
            else: raise ValueError

            for s in srange:
                for a in arange:
                    lo.append(lo_); hi.append(hi_)
                    coefs.append(flatten([[LegendreProduct(s,a,i,j)
                                            for j in range(self.size)]
                                                for i in range(self.size)]))
        return coefs, lo, hi

    def submit(self, root : str) -> None:
        '''Create a directory for the fit and populate files needed to run'''

        pth = join(root, self.name)
        safeMkdir(pth)

        if not exists(join(pth, 'result.json')):
            train = join(pth, 'train.py')
            copy('/Users/ksb/functionals/functionals/scripts/train.py', train)
            files = ['params',  'data',     'constraint']
            data  = [self.args, self.data(), self.const()]
            for fname,d in zip(files, data):
                with open(join(pth, fname+'.json'), 'w') as f:
                    dump(d, f)
        else:
            print('Already done')
