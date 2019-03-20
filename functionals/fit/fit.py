# External
from typing     import Dict as D, List as L, Tuple as T, Any, Optional as O
from ast        import literal_eval
from time       import time
from os.path    import join
from shutil     import copyfile
from os         import system
from json       import dumps, loads, dump, load
from io         import StringIO
from contextlib import redirect_stdout
from traceback  import format_exc
from contextlib import contextmanager
from hashlib import sha512
from numpy.linalg              import lstsq,norm                            # type: ignore

import sys

from scipy.optimize import least_squares # type: ignore

from numpy          import (average,empty,tile,eye,sum,ones,array,linspace,inf,  # type: ignore
                            zeros,vstack,concatenate as concat)
from numpy.linalg   import lstsq,inv                            # type: ignore
import warnings; warnings.filterwarnings("ignore")
from scipy.stats import linregress         # type: ignore
from MySQLdb     import connect,Connection # type: ignore

# Internal
from functionals.fit.utilities import LegendreProduct
from functionals.fit.data import process_data,weight
from functionals.scripts.fit.h_norm_const import intkernel,quad

################################################################################
# Helper functions
#-----------------
def flatten(lol:list)->list: return [item for sublist in lol for item in sublist]

def sqlselect(conn : Connection, q : str, binds : list = []) -> L[tuple]:
    with conn.cursor() as cxn: # type: ignore
        cxn.execute(q,args=binds)
        return cxn.fetchall()

def hash_(x:Any)->str:
    return sha512(str(x).encode()).hexdigest()[:10]

def h_norm_vec(a1:float,msb:float)->array:
    return [quad(intkernel(i,j,a1,msb), 0., inf)[0] for j in range(8) for i in range(8)]


# SQL
q1 = '''SELECT const_name FROM const'''

q2 = '''SELECT dataconst,bm_weight,lat_weight,reg,constden,consts
        FROM fitparams WHERE fitparams_id=%s'''

q3 = '''SELECT expt.n_atoms, species.n_atoms, volumes, energies,
                contribs, expt.composition, atomic_contribs,
                atomic_energies, bulk_contribs, bulk_energy,
                expt_cohesive_energy, expt_bm, expt_volume
        FROM expt JOIN species ON species=species_id
        WHERE calc = %s
              AND name REGEXP BINARY(%s)'''

q4 = '''SELECT pw,econv,data,a11,a12,a13,a14,a15,msb
        FROM calc JOIN functional ON functional=functional_id
        WHERE calc_id = %s'''

q5 = '''SELECT const_name,val,kind,s,alpha
        FROM const'''

q7 = '''SELECT expt.n_atoms,expt.composition,symmetry
        FROM expt JOIN species ON species=species_id
        WHERE calc = %s AND name REGEXP BINARY(%s)'''

class Fit(object):
    """
    Fit a 8 x 8 matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """

    def __init__(self,
                 constraints    : L[dict],
                 calc           : dict, # calc id
                 decay          : int,  # which a1 output are we using for the data
                 bmw            : float,
                 lw             : float,
                 reg            : float,
                 cd             : float,
                 a1             : float,
                 msb            : float,
                 dataconstraint : str,
                 ceX            : array,
                 bmX            : array,
                 lX             : array,
                 ceY            : array,
                 bmY            : array,
                 lY             : array,
                 expts : L[T[str,str,str]],
                 ) -> None:
        self.constraints = constraints
        self.calc        = calc

        self.dataconstraint = dataconstraint
        self.decay       = decay
        self.bmw         = bmw
        self.lw          = lw
        self.lam_reg     = reg
        self.cd          = cd
        self.a1         = a1
        self.msb = msb
        self.ceX,self.bmX,self.lX,self.ceY,self.bmY,self.lY=ceX,bmX,lX,ceY,bmY,lY
        self.data = [ceX,bmX,lX,ceY,bmY,lY]
        self.X,self.Y = vstack((ceX,bmX,lX)), concat((ceY,bmY,lY))
        self.expts = expts

        self.c_A_eq,self.c_b_eq,self.c_A_lt,self.c_b_lt = self.const_xy()
        assert decay in range(5)

        self.kwargs = dict(method        = 'trf',
                           ftol          = 1e-12,
                           xtol          = 1e-12,
                           gtol          = 1e-12,
                           x_scale       = 'jac',
                           loss          = 'linear',
                           max_nfev      = 100000,
                           verbose       = 1) # type: dict
    ###############
    # MAIN METHOD #
    ###############
    def fit(self)->T[float,list]:
        '''Returns: time in seconds to fit, trajectory of solutions until convergence'''
        # Start fitting
        ################

        # initialize
        start_time = time()
        x0         = lstsq(self.X,self.Y,rcond=None)[0] #
        bounds     = ([min(x0)] * 64, [max(x0)]* 64)
        dx         = float('inf')
        c_reg      = self.lam_reg/10000 # initial regularization penalty = 1/1000 of final
        constr     = 1e-10               # negligible constraints initially
        sols       = []

        # main loop - break when constraints are high AND negligible change in x
        while constr < 1e5 or dx > self.kwargs['xtol']:
            # Do fitting
            kwargs = dict(constr=constr,c_reg=c_reg)
            res = least_squares(self.cost, x0, jac = self.dcost, bounds=bounds, kwargs = kwargs,**self.kwargs)

            # Do updates
            sols.append(x0.tolist())
            dx      = norm(res.x - x0)            # how much has solution changed
            x0      = res.x                       # overwrite previous solution
            c_reg  += (self.lam_reg - c_reg) / 5 # asymptotically increase regularization
            constr *= 5                           # increase constraints

        return time() - start_time, sols

    @classmethod
    def costs(cls, db : str, x_ : str, calc : int, decay : int) -> T[float,float,float,float]:
        from numpy.linalg import norm

        def r2(vec:array,A:array,b:array)->float:
            yhat = A @ vec
            _,_,r,_,_  = linregress(yhat,b)
            return r**2

        x = loads(x_)
        fit = cls.from_db(db,fp_id=None,calc_id=calc,decay=decay)
        ceX,bmX,lX,ceY,bmY,lY = fit.data

        ce,bm,l = [r2(x,array(getattr(fit,a+'X')),array(getattr(fit,a+'Y')))
                    for a in ['ce','bm','l']]
        c_viol = norm(fit.c_loss(x))

        return ce,bm,l,c_viol

    ################
    # CONSTRUCTORS #
    ################
    @classmethod
    def from_db(cls,
                db          : str,
                fp_id       : O[int],
                calc_id     : int,
                decay       : int,
               )->'Fit':
        with open(db,'r') as fi: conn = connect(**load(fi),  autocommit  = True)

        # Extract fitting parameters
        if fp_id is not None:
            dc,bmw_,lw_,reg_,cd,constrs_ = sqlselect(conn,q2,[fp_id])[0]
            constraints = constrs_.split()
            bmw,lw,reg = map(float,[bmw_,lw_,reg_])
        else:

            bmw,lw,reg = 1,1,0 #these don't matter
            constraints = [x for (x,) in sqlselect(conn,q1)] # WE ARE DOING A GLOBAL ANALYSIS: INCLUDE EVERY CONSTRAINT
            dc,cd = '.',10 # include all data, fine grid density
        # Extract constraints
        con = []
        for name,val,kind,s_,alpha in sqlselect(conn,q5):
            if name in constraints:
                s,a,v   = [float(x) if x is not None else None for x in
                                [s_,alpha,val]]
                w = len(constraints) - constraints.index(name)
                con.append(dict(name=name,weight=w,s=s,alpha=a,kind=kind,val=v))

        pw_,econv_,fx_,a1_,a2_,a3_,a4_,a5_,msb_ = sqlselect(conn,q4,[calc_id])[0]
        pw,econv,a11,a12,a13,a14,a15,msb = map(float,[pw_,econv_,a1_,a2_,a3_,a4_,a5_,msb_])
        calc = dict(pw=pw,econv=econv,a11=a11,a12=a12,a13=a13,a14=a14,a15=a15,msb=msb,fx=loads(fx_))
        a1 = [a11,a12,a13,a14,a15][decay]
        expts = cls.get_expts(conn,calc_id,dc)
        ceX,bmX,lX,ceY,bmY,lY = cls.data_xy(db=conn,calc=calc_id,fx=calc['fx'],decay=decay,dconstraint=dc,bmw=bmw,lw=lw)
        return cls(constraints=con,calc=calc,decay=decay,cd=cd,bmw=bmw,lw=lw,dataconstraint=dc,
                   reg=reg,a1=a1,msb=msb,ceX=ceX,bmX=bmX,lX=lX,ceY=ceY,bmY=bmY,lY=lY,expts=expts)

    @classmethod
    def from_json(cls,pth : str)->'Fit':
        with open(join(pth,'constraints.json'),'r') as fi: con = load(fi)
        with open(join(pth,'metadata.json'),'r')    as fi: md = load(fi)
        with open(join(pth,'data.json'),'r')        as fi: d = load(fi)
        p = md['params']
        ceX,bmX,lX,ceY,bmY,lY = map(array,d)
        return cls(constraints=con,calc=md['calc'],decay=p['decay'],cd=p['cd'],
                    bmw=p['bmw'],lw=p['lw'],reg=p['reg'],a1=p['a1'],msb=p['msb'],
                    dataconstraint=p['dc'],ceX=ceX,bmX=bmX,lX=lX,ceY=ceY,bmY=bmY,lY=lY,
                    expts=md['expts'])

    #######################
    # COST/LOSS FUNCTIONS #
    #######################
    # Constraints
    @staticmethod
    def relu(x:array) -> array: return x * (x > 0)
    def c_lt_loss(self,x : array) -> array:  return self.relu(self.c_A_lt @ x - self.c_b_lt)
    def jc_lt_loss(self,x : array) -> array: return self.c_A_lt * tile(((self.c_A_lt @ x - self.c_b_lt) > 0).reshape(-1,1),(1,64))

    def c_eq_loss(self,x  : array) -> array: return self.c_A_eq @ x - self.c_b_eq
    def jc_eq_loss(self,_ : array) -> array: return self.c_A_eq

    def c_loss(self,x  : array) -> float: return concat((self.c_lt_loss(x),  self.c_eq_loss(x)))
    def dc_loss(self,x : array) -> array: return vstack((self.jc_lt_loss(x),self.jc_eq_loss(x)))

    # Regularization
    @staticmethod
    def reg(x  : array) -> array: return x
    @staticmethod
    def dreg(_ : array) -> array: return eye(64)

    # Data
    def loss(self,x  : array) -> array: return self.X @ x - self.Y
    def dloss(self,_ : array) -> array: return self.X

    # All together: data, regularization, constraints
    def cost(self,x:array,constr:float,c_reg:float)->array:
        return concat(( self.loss(x),constr * self.c_loss(x),c_reg * self.reg(x)))

    def dcost(self,x:array,constr:float,c_reg:float)->array:
        return vstack((self.dloss(x),constr * self.dc_loss(x),c_reg * self.dreg(x)))

    ##########
    # HELPER #
    ##########
    def _linconst(self,const : dict) -> T[bool,array,array]:
        ''' Turns a linear constraint into an Ax = b matrix '''
        w = const['weight']
        assert w > 0
        grid = linspace(0,5,self.cd).tolist()

        val = float(const['val'])

        if   const['kind'] == 'eq':  eq = True  ; sign = 1
        elif const['kind'] == 'lt':  eq = False ; sign = 1
        elif const['kind'] == 'gt':  eq = False ; sign = -1
        else: raise ValueError

        s_, a_ = const['s'],const['alpha']
        srange = [float(s_)] if s_ is not None else grid
        arange = [float(a_)] if a_ is not None else grid
        card   = len(srange) * len(arange)
        w      = w / card # reduce weight if many data points
        vals   = [val for _ in range(card)]
        c_A    = [flatten([[LegendreProduct(s,a,self.a1,self.msb,i,j)
                            for j in range(8)] for i in range(8)])
                        for s in srange for a in arange]
        return eq, w * array(c_A)  * sign, w * array(vals) * sign


    def const_xy(self)->T[array,array,array,array]:
        '''
        Convert constraint dictionaries into matrices (add analytic H solution, too)
        '''
        c_A_eq,c_b_eq,c_A_lt,c_b_lt = empty((0,64)),empty(0),empty((0,64)),empty(0)
        for const in self.constraints:
                eq,A,b = self._linconst(const)
                if eq:
                    c_A_eq = vstack((c_A_eq,A))
                    c_b_eq = concat((c_b_eq,b))
                else:
                    c_A_lt = vstack((c_A_lt,A))
                    c_b_lt = concat((c_b_lt,b))

        c_A_eq = vstack((c_A_eq,h_norm_vec(self.a1,self.msb)))
        c_b_eq = concat((c_b_eq,[-.3125]))
        return c_A_eq,c_b_eq,c_A_lt,c_b_lt

    def metadata(self) -> dict:
        '''Summary of parameters used to determine fit inputs'''
        c = [x['name'] for x in self.constraints]
        return dict(expts = self.expts, cons=[c['name'] for c in self.constraints],
                    params = dict(bmw=self.bmw,lw=self.lw,cd=self.cd,c=c,decay=self.decay,
                                reg=self.lam_reg,dc=self.dataconstraint,a1=self.a1,msb=self.msb),
                    calc   = self.calc)

    @staticmethod
    def get_expts(db:str,calc:int,dataconst:str) -> L[T[str,str,str]]:
        '''Identifying info about which DFT jobs were used in the fit'''
        return sqlselect(db,q7,[calc,dataconst]) # type: ignore

    @classmethod
    def data_xy(cls,
                db          : Connection,
                calc        : int,
                fx          : array,
                decay       : int,
                dconstraint : str,
                bmw         : float,
                lw          : float
               ) -> T[array,array,array,array,array,array]:
            sli = slice(64*decay,64*(decay+1))
            fxi = fx[sli]
            ceX_,bmX_,lX_,ceY_,bmY_,lY_ = [],[],[],[],[],[]
            for n,nsp,vols,engs,xcs,comp,ac,ae,bc,be,ce,bm,vol in sqlselect(db,q3,[calc,dconstraint]):
                a_contribs = {k:array(v)[sli] for k,v in literal_eval(ac).items()}
                b_contribs = array(loads(bc)[sli])
                contribvec = [xc[sli] for xc in loads(xcs)]
                new_cex,new_cey = cls.cohesive(comp,a_contribs,ae,fx,n//nsp,ce,be,b_contribs)
                ceX_.append(new_cex);ceY_.append(new_cey)
                new_bmx,new_bmy,new_lx,new_ly = cls.bm_lat(engs,vols,contribvec,float(bm),float(vol),fx)
                lX_.append(new_lx);bmX_.append(new_bmx)
                bmY_.append(new_ly);lY_.append(new_bmy)

            ceX,bmX,lX = [array(x) if x else empty((0,64)) for x in [ceX_,bmX_,lX_]]
            ceY,bmY,lY = [array(x) if x else empty(0) for x in [ceY_,bmY_,lY_]]

            ce_avg,bm_avg,l_avg = [average(x) if 0 not in x.shape else 0 for x in [ceY,bmY,lY]]
            pre_l = lw * ce_avg / l_avg
            pre_b = bmw * ce_avg / bm_avg
            return (ceX, pre_b * bmX, pre_l * lX,
                    ceY, pre_b * bmY, pre_l * lY)

    def write(self, pth : str)->None:
        root = '/Users/ksb/functionals/functionals/scripts/fit/'
        def write(fi:str,x:Any)->None:
            with open(join(pth,fi)+'.json','w') as file: dump(x,file)

        # Write to directory
        #--------------------
        md = self.metadata()
        md['uid'] = self.uid()
        write('metadata', md)
        write('constraints', self.constraints)
        write('data', [x.tolist() for x in self.data])
        copyfile(root+'runfit.py',join(pth,'runfit.py'))
        copyfile(root+'subfit.sh',join(pth,'subfit.sh'))
        system('chmod 755 '+join(pth,'subfit.sh'))

    def uid(self) -> str:
        md = self.metadata()
        md['params'].pop('decay'); md['params'].pop('a1') # we want things w/ different decay to have same 'name'
        return hash_(md)

    @staticmethod
    def cohesive(comp_   : str,
                 raw_ac_ : D[int,array],
                 raw_ae_ : str,
                 coefs   : array,
                 ratio   : int,
                 tar     : float,
                 ebulk_  : float,
                 bcontribs_: array,
                ) -> T[list,list]:
        '''
        Convert a row data into a row of coefficients and a target for cohesive data
        fitting. Requires data to have the following keys:
        - composition, atomic_contribs, atomic_energies, coefs, bulk_ratio,
          bulk_energy, bulk_contribs, expt_cohesive_energy

        returns a 5x64 matrix and a length-5 vector
        '''

        # Extract info from dictionary
        comp   = literal_eval(comp_)             # type: D[int,int]
        raw_ac = raw_ac_.items() # int -> 64 element array
        raw_ae = literal_eval(raw_ae_).items() # int -> float

        ex_ce  = float(tar)           # experimental E atom - E bulk
        e_bulk = float(ebulk_) / ratio            # calculated energy of bulk reference system

        contribs = bcontribs_.flatten()

        x_bulk = contribs / ratio # 64 element array

        # Analysis
        #---------
        # Get coefficients representing the change in exchange contributions
        atom_contribs = {i:array(xs) for i,xs in raw_ac} # int -> 8x8 matrix
        atom_energies = {i:float(e)  for i,e  in raw_ae} # int -> float

        e_atom  = sum([atom_energies[e]*num for e,num in comp.items()]) # float
        x_atom  = sum([atom_contribs[e]*num for e,num in comp.items()],axis=0) # 8x8 matrix

        dx = (x_atom - x_bulk) # 64 element array, THIS IS WHAT WE ARE FITTING

        # Get target to fit the above to: JUST the ex_component of cohesive energy
        ex_atom = coefs @ x_atom # Use the BEEF coefficients
        ex_bulk = coefs @ x_bulk # from this particular calculator

        nonx_e_atom = e_atom - ex_atom
        nonx_e_bulk = e_bulk - ex_bulk
        target      = ex_ce  - (nonx_e_atom - nonx_e_bulk) # Just Ex_atom - Ex_bulk


        return dx.tolist(),target

    @staticmethod
    def bm_lat(engs     : str,
               vols     : str,
               contribs : list,
               expt_bm  : float,
               expt_vol : float,
               coefs    : array,
              ) -> T[array,array,array,array]:
        '''
        Fit energies to quadratic form using linear algebra
        ---------------------------------------------------
        Solve b = A x
            with: x = (Aᵗ·A)⁻¹ · Aᵗ · b
                - x = vector of three elements: x2,x1,x0
                    - such that: Energy = x2*Vol²+x1*Vol+x0
                - b = contribs·coefs + e_nonx
                    - this is a function of the vector being fit
                - A = a "vandermonde" matrix generated purely by the volume data

        so x = (Aᵗ·A)⁻¹·Aᵗ·(contribs·coefs + e_nonx)
            or x[i](coefs) = VECTOR · COEFS + CONST
                - VECTOR = i'th row of (Aᵗ·A)⁻¹·Aᵗ·contribs
                - CONST  = i'th element of 1 x 3 vector: (Aᵗ·A)⁻¹·Aᵗ·contribs · e_nonx

        Curvature prediction error = (x[2](coefs) - expt_curv)
        Lattice prediction error   = (x[1](coefs) / 2*expt_curv) - expt_volume
            - note that minimum of parabola is at -x1/2*x2

        In order to make the error some A·x - b, we need:
            Curvature:
                - A = x[2] VECTOR
                - b = expt_curv - x[2] CONST
            Lattice:
                - A = x[1] VECTOR / 2*expt_curv
                - b = expt_vol  - (x[1] CONST / 2*expt_curv)

        We need this for all 5 functionals, so our return types are
        A_bm --- 5 x 64
        '''

        # Common stuff to preprocessing both BM and Lattice data
        #----------------------------------------------------------

        energies = array(loads(engs))   # 5 element array
        volumes  = array(loads(vols))   # 5 element array

        contribs = array(contribs).reshape((5,64)) # 5 x 64

        bm          = expt_bm*10**9     # experimental value, Pa or N / m²
        expt_volume = expt_vol*(10**-30)       # experimental volume, m^3

        # Get experimental curvature to E vs V
        #--------------------------------------
        curv_      = bm / expt_volume                   # experimental d²E/dV², J/m^6 = N/m^5
        expt_curv  = curv_ * (10**-60) * (6.242*10**18) # eV / A^6


        e_nonx  = energies - contribs @ coefs # len-5 vector

        # also because we only need "a", we only dot the last row with the coef (col) vec

        vander = vstack((ones(len(volumes)),volumes,volumes**2)).T  # 5 x 3
        vinv   = inv(vander.T @ vander)                             # 3 x 3
        solver =  vinv @ vander.T                                   # 3 x 5
        vecs    = solver @ contribs                                        # 3 x 64*5
        constvec  = solver @ e_nonx                                 # 1 x 3

        curv_vec   = vecs[2]
        curv_const = expt_curv - constvec[2]
        lat_vec    = vecs[1] / 2*expt_curv
        lat_const  = expt_volume -  constvec[1] / 2*expt_curv

        return  curv_vec.tolist(), curv_const, lat_vec.tolist(), lat_const


"""
BONUS
-----

Derivation of Jacobian in Einstein notation:

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

Derivation of constant hessian:
        To see this, take derivative of jac result w/r/t some Rm:

        d(FjFj)/d(Rk)d(Rm) = 2 * [A]jk * d(Fj)/d(Rm) =  2 * [A]jk * [A]jm
"""
