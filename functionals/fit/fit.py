# External
from typing       import Dict as D, List as L, Tuple as T, Any, Optional as O, Callable as C
from ast          import literal_eval
from time         import time
from os.path      import join
from copy         import deepcopy
from sys          import stdout
from shutil       import copyfile
from os           import system, environ
from json         import dumps, loads, dump, load,JSONEncoder
from io           import StringIO
from contextlib   import redirect_stdout
from traceback    import format_exc
from contextlib   import contextmanager
from hashlib      import sha512

from numpy          import (ndarray,average,empty,tile,eye,sum,ones,array,linspace,inf,  # type: ignore
                            zeros,inf,vstack,concatenate as concat,diag, multiply,
                            stack, array_equal)
from numpy.linalg   import lstsq,norm,inv,LinAlgError      # type: ignore

from scipy.optimize import least_squares # type: ignore

#import warnings; warnings.filterwarnings("ignore")

from scipy.stats import linregress         # type: ignore
from psycopg2    import connect # type: ignore

# Internal
from functionals.fit.constraint           import Constraint,lda,liebox,scan11,pos,hnorm
from functionals.fit.data                 import process_data,weight

Connection = Any
Arrs       = L[array]
################################################################################
# Helper functions
#-----------------

def sqlselect(conn : Connection, q : str, binds : list = []) -> L[tuple]:
    with conn.cursor() as cxn: # type: ignore
        cxn.execute(q,vars=binds)
        return cxn.fetchall()

def hash_(x : Any)->str: return sha512(str(x).encode()).hexdigest()[:10]


bias = array([10**(x%8 + x//8) for x in range(64)]) # unequally weigh matrix elements

def safeLoad(x:Any)->Any: return x if x is None else loads(x)

q2 = '''SELECT bm_weight,lat_weight,reg,consts
        FROM fitparams WHERE fitparams_id=%s'''

q3 = '''SELECT a_ce,a_bm,a_l,b_ce,b_bm,b_l,expt_ce,expt_bm,expt_vol,name
        FROM bulks JOIN job on job=job_id
        WHERE calc = %s'''

q4 = '''SELECT pw,econv,data,a11,a12,a13,a14,a15,msb
        FROM calc C JOIN beef B on B.functional = C.functional
        WHERE calc_id = %s'''

class Fit(object):
    """
    Fit a 8 x 8 matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """

    def __init__(self,
                 cons : L[str],
                 calc : dict,
                 bmw  : float, lw   : float, reg  : float,
                 ceX : list, bmX : list, lX : list,
                 ceB : list, bmB : list, lB : list,
                 ceY : list, bmY : list, lY : list,
                 pre_ce : float, pre_bm : float, pre_l : float,
                 c_A_eq: list,c_b_eq: list, c_A_lt: list, c_b_lt: list,
                 C_A_eq: list,C_b_eq: list, C_A_lt: list, C_b_lt: list,
                 X:list,Y:list) -> None:

        self.cons = cons
        self.calc = calc
        self.bmw = bmw; self.lw = lw; self.lam_reg = reg;

        self.pre_ce,self.pre_bm,self.pre_l = pre_ce, pre_bm, pre_l


        self.ceY, self.bmY, self.lY,\
        self.c_b_eq,self.c_b_lt,self.C_b_eq,self.C_b_lt\
            = map(array, [ceY, bmY, lY,c_b_eq,c_b_lt,C_b_eq,C_b_lt])

        self.ceX, self.bmX, self.lX, self.ceB, self.bmB, self.lB, \
        self.X,self.Y \
          = [[array(x[i]) for i in range(5)] for x in
            [ceX, bmX, lX, ceB, bmB, lB,
             X,Y]]

        self.c_A_eq, self.c_A_lt, self.C_A_eq, self.C_A_lt,\
          = [[array(x[i]) if x[i] else empty((0,64)) for i in range(5)] for x in
            [c_A_eq, c_A_lt, C_A_eq, C_A_lt]]

        self.kwargs = dict(method   = 'trf',
                           ftol = 1e-12, xtol = 1e-12, gtol = 1e-12,
                           x_scale  = 'jac',
                           loss     = 'linear',
                           max_nfev = 50000,
                           verbose  = 0) # type: dict


    def __eq__(self, other : object) -> bool:
        return False if not isinstance(other,Fit) else \
             all([self.calc         == other.calc,
                  self.cons         == other.cons,
                  self.bmw          == other.bmw,
                  self.lw           == other.lw,
                  self.lam_reg      == other.lam_reg,
                  self.pre_ce       == other.pre_ce,
                  self.pre_bm       == other.pre_bm,
                  self.pre_l        == other.pre_l,
                  array_equal(self.ceY, other.ceY),
                  array_equal(self.bmY, other.bmY),
                  array_equal(self.lY, other.lY),
                  array_equal(self.ceX,other.ceX),
                  array_equal(self.bmX,other.bmX),
                  array_equal(self.lX,other.lX),
                  array_equal(self.ceB,other.ceB),
                  array_equal(self.bmB,other.bmB),
                  array_equal(self.lB,other.lB),
                  array_equal(self.c_A_eq,other.c_A_eq),
                  array_equal(self.c_A_lt,other.c_A_lt),
                  array_equal(self.C_A_eq,other.C_A_eq),
                  array_equal(self.C_A_lt,other.C_A_lt),
                  array_equal(self.Y,other.Y),
                  array_equal(self.X,other.X)])

    @property
    def beef(self)->array:
        return array([
            1.18029330e+00,8.53027860e-03,-1.02312143e-01,6.85757490e-02,-6.61294786e-03,-2.84176163e-02,5.54283363e-03,3.95434277e-03,
            -1.98479086e-03,1.00339208e-01,-4.34643460e-02,-1.82177954e-02,1.62638575e-02,-8.84148272e-03,-9.57417512e-03,9.40675747e-03,
            6.37590839e-03,-8.79090772e-03,-1.50103636e-02,2.80678872e-02,-1.82911291e-02,-1.88495102e-02,1.69805915e-07,-2.76524680e-07,
            1.44642135e-03,-3.03347141e-03,2.93253041e-03,-8.45508103e-03,6.31891628e-03,-8.96771404e-03,-2.65114646e-08,5.05920757e-08,
            6.65511484e-04,1.19130546e-03,1.82906057e-03,3.39308972e-03,-7.90811707e-08,1.62238741e-07,-4.16393106e-08,5.54588743e-08,
            -1.16063796e-04,8.22139896e-04,-3.51041030e-04,8.96739466e-04,2.09603871e-08,-3.76702959e-08,2.36391411e-08,-3.38128188e-08,
            -5.54173599e-06,-5.14204676e-05,6.68980219e-09,-2.16860568e-08,9.12223751e-09,-1.38472194e-08,6.94482484e-09,-7.74224962e-09,
            7.36062570e-07,-9.40351563e-06,-2.23014657e-09,6.74910119e-09,-4.93824365e-09,8.50272392e-09,-6.91592964e-09,8.88525527e-09])

    all_cons = dict(lda    = lda,
                    liebox = liebox,
                    scan11 = scan11,
                    pos    = pos,
                    hnorm  = hnorm) # type: D[str,Constraint]

    @property
    def constraints(self)->L[Constraint]:
        return [self.all_cons[x] for x in self.cons]

    ###############
    # MAIN METHOD #
    ###############
    def fit(self)->list:
        '''Returns: trajectory of solutions until convergence'''

        # initialize
        start_time = time()
        allsols    = []

        for i in range(5):
            #print('\n\ti = %d'%i,end='')
            x0         = lstsq(self.X[i],self.Y[i],rcond=None)[0] #array([1]+[0]*63) #
            dx         = float('inf')
            c_reg      = self.lam_reg/10000 # initial regularization penalty = 1/1000 of final
            constr     = 1e-10               # negligible constraints initially
            sols       = [[self.beef.tolist(), self.r2avg(self.beef,i),norm(self.c_loss(self.beef,i))]]
            bounds     = ([-dx]*64,[dx]*64)#([-2.]+[-0.5] * 63, [2.]+[0.5]* 63)

            # main loop - break when constraints are high AND negligible change in x
            while (constr < 1e5) or (dx > self.kwargs['xtol']):
                # Do fitting
                kwargs = dict(constr=constr,c_reg=c_reg)
                res = least_squares(self.cost(i), x0, jac = self.dcost(i), bounds=bounds, kwargs = kwargs,**self.kwargs)#try: except LinAlgError:
                # Do updates
                sols.append([x0.tolist(),self.r2avg(x0,i),norm(self.c_loss(x0,i))])
                dx      = norm(res.x - x0)            # how much has solution changed
                x0      = res.x                       # overwrite previous solution
                c_reg  += (self.lam_reg - c_reg) / 5 # asymptotically increase regularization
                constr *= 1.01                           # increase constraints
            allsols.append(sols)
        tottime = time() - start_time # currently do nothing with this
        assert len(allsols)==5
        return allsols

    def costs(self, x : Arrs) -> T[L[float],L[float],L[float],L[float],L[float],L[float]]:
        '''Returns the following metrics: MSE - cohesive (eV), bulkmod (GPa), lattice (A), three r2 values, and constraint violation'''
        from numpy.linalg import norm

        yys = [([ce@X+b for X,b,ce in zip(x,self.ceB,self.ceX)],   self.ceY),
               ([bm@X+b for X,b,bm in zip(x,self.bmB,self.bmX)],   self.bmY),
               ([lx@X+b for X,b,lx in zip(x,self.lB,self.lX)],     self.lY)] # convert volume -> lattice constant

        (rces,ces),(rbms,bms),(rls,ls) = [([linregress(y_,y)[2]**2 for y_ in yy_],
                                           [average((y_ - y)**2)    for y_ in yy_])
                                            for yy_,y in yys]

        return ces,bms,ls,rces,rbms,rls


    def midcosts(self, x : Arrs) -> T[float,float,float,float,float,float]:
        '''Returns the cost metrics but only for the central functional (a1 = a13)'''
        ce,bm,l,rce,rbm,rl = [x[2] for x in self.costs(x)]
        return ce,bm,l,rce,rbm,rl

    def decaycosts(self,  x : str) -> T[float,float,float,float,float]:
        '''Avg R2 value for the 5 decays'''
        x_ = [array(y) for y in loads(x)]
        _,_,_,rce,rbm,rl = self.costs(x_)
        x1,x2,x3,x4,x5 = [(r1+r2+r3)/3 for r1,r2,r3 in zip(rce,rbm,rl)]
        return x1,x2,x3,x4,x5

    def allcosts(self, x : str) -> T[str,float,float,float,float,float,float]:
        '''Combine midcosts and decay costs information'''
        x_ = [array(y) for y in loads(x)]
        ces,bms,ls,rces,rbms,rls = self.costs(x_)
        x1,x2,x3,x4,x5 = [(r1+r2+r3)/3 for r1,r2,r3 in zip(rces,rbms,rls)]
        ce,bm,l,rce,rbm,rl = [min(x[2],10**13) for x in [ces,bms,ls,rces,rbms,rls]]
        return dumps([x1,x2,x3,x4,x5]),ce,bm,l,rce,rbm,rl

    def allviols(self, x : array, decay : int) -> D[str,float]:
        '''Report all constraint violation magnitudes'''
        out = dict(tot=0.)
        for k,v in self.all_cons.items():
            out[k]  = v.viol(x,self.calc['decays'][decay],self.calc['msb'])
            if k in self.cons:
                out['tot']       += out[k] * v.weight
        return out

    ################
    # CONSTRUCTORS #
    ################
    @classmethod
    def from_db(cls,
                db      : str,
                fp_id   : int,
                calc_id : int,
               )->'Fit':
        with open(db,'r') as fi:
            kwargs = load(fi)
            kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
            conn = connect(**kwargs)

        # Extract fitting parameters
        bmw_,lw_,reg_,constrs_ = sqlselect(conn,q2,[fp_id])[0]
        constraints = constrs_.split()
        bmw,lw,reg = map(float,[bmw_,lw_,reg_])

        # Extract calculation parameters
        pw_,econv_,fx_,a1_,a2_,a3_,a4_,a5_,msb_ = sqlselect(conn,q4,[calc_id])[0]
        pw,econv,a11,a12,a13,a14,a15,msb = map(float,[pw_,econv_,a1_,a2_,a3_,a4_,a5_,msb_])
        calc = dict(pw=pw,econv=econv,decays=[a11,a12,a13,a14,a15],msb=msb,fx=loads(fx_))

        # Extract fitted data
        ceX,bmX,lX,ceB,bmB,lB,ceY,bmY,lY,pre_ce,pre_bm,pre_l = cls.data_xy(db=conn,calc=calc_id,fx=calc['fx'],bmw=bmw,lw=lw)

        c_A_eq,c_b_eq, c_A_lt,c_b_lt = Constraint.const_xy([cls.all_cons[x] for x in constraints],msb,calc['decays'])
        C_A_eq,C_b_eq,C_A_lt,C_b_lt = Constraint.const_xy(list(cls.all_cons.values()),msb,calc['decays'])

        X = [vstack((pre_ce * array(xc),
                     pre_bm * array(xb),
                     pre_l  * array(xl))).tolist()
             for xc,xb,xl in  zip(ceX,bmX,lX)]

        Y = [concat((pre_ce * (array(ceY) - array(bc)),
                     pre_bm * (array(bmY) - array(bb)),
                     pre_l  * (array(lY)  - array(bl)))).tolist()
             for bc,bb,bl in zip(ceB,bmB,lB)]

        return cls(cons= constraints, calc=calc, bmw=bmw, lw=lw,
                   reg=reg, ceX=ceX, bmX=bmX, lX=lX,ceB=ceB, bmB=bmB, lB=lB,
                   ceY=ceY,bmY=bmY,lY=lY,
                   pre_ce=pre_ce, pre_bm=pre_bm, pre_l=pre_l,
                   c_A_eq=c_A_eq, c_b_eq=c_b_eq, c_A_lt=c_A_lt,c_b_lt=c_b_lt,
                   C_A_eq=C_A_eq, C_b_eq=C_b_eq, C_A_lt=C_A_lt,C_b_lt=C_b_lt,
                   X=X,Y=Y)

    @classmethod
    def from_json(cls,pth : str)->'Fit':
        with open(join(pth,'metadata.json'),'r')    as fi: md = load(fi)
        with open(join(pth,'data.json'),'r')        as fi: d = load(fi)
        p = md['params']
        ceX,bmX,lX,ceB,bmB,lB,ceY,bmY,lY,c_A_eq,c_b_eq, c_A_lt,c_b_lt,C_A_eq,C_b_eq,C_A_lt,C_b_lt,X,Y = d

        return cls(cons=p['c'],calc=md['calc'],
                    bmw=p['bmw'],lw=p['lw'],reg=p['reg'],
                    ceX=ceX,bmX=bmX,lX=lX,ceY=ceY,bmY=bmY,lY=lY,ceB=ceB,bmB=bmB,lB=lB,
                    c_A_eq=c_A_eq,c_b_eq=c_b_eq, c_A_lt=c_A_lt,c_b_lt=c_b_lt,
                    C_A_eq=C_A_eq,C_b_eq=C_b_eq, C_A_lt=C_A_lt,C_b_lt=C_b_lt,
                    pre_ce=p['pre_ce'],pre_bm=p['pre_bm'],pre_l=p['pre_l'],X=X,Y=Y)

    #######################
    # COST/LOSS FUNCTIONS #
    #######################
    # Constraints
    @staticmethod
    def relu(x:array) -> array: return x * (x > 0)
    def c_lt_loss(self,x : array, i : int, all : bool = False) -> array:
        if all: return self.relu(self.C_A_lt[i] @ x - self.C_b_lt)
        else:   return self.relu(self.c_A_lt[i] @ x - self.c_b_lt)
    def jc_lt_loss(self,x : array, i : int) -> array: return self.c_A_lt[i] * tile(((self.c_A_lt[i] @ x - self.c_b_lt) > 0).reshape(-1,1),(1,64))

    def c_eq_loss(self,x  : array, i : int, all : bool = False) -> array:
        if all: return self.C_A_eq[i] @ x - self.C_b_eq #use all known constraints
        else:   return self.c_A_eq[i] @ x - self.c_b_eq
    def jc_eq_loss(self,_ : array, i : int) -> array: return self.c_A_eq[i]

    def c_loss(self,x  : array, i : int, all : bool = False) -> array:
        return concat((self.c_lt_loss(x,i,all),  self.c_eq_loss(x,i,all)))
    def dc_loss(self,x : array, i : int) -> array:
        return vstack((self.jc_lt_loss(x,i),self.jc_eq_loss(x,i)))

    # Regularization
    @staticmethod
    def reg(x  : array, _ : int) -> array: return multiply(bias,x)
    @staticmethod
    def dreg(_ : array, __ : int) -> array: return diag(bias)

    # Data
    def loss(self,x  : array, i : int) -> array: return self.X[i] @ x - self.Y[i]
    def dloss(self,_ : array, i : int) -> array: return self.X[i]

    # All together: data, regularization, constraints
    def cost(self, i : int) -> C:
        def f(x : array, constr : float, c_reg : float) -> array:
            return concat((self.loss(x, i), constr * self.c_loss(x, i), c_reg * self.reg(x, i)))
        return f

    def dcost(self, i : int) -> C :
        def f(x : array, constr : float, c_reg : float) -> array:
            return vstack((self.dloss(x, i), constr * self.dc_loss(x, i), c_reg * self.dreg(x, i)))
        return f

    # Data
    def r2avg(self,x  : array, i : int) -> array:
        '''Weighted average R2 value for an arbitrary BEEF vector'''
        ys = [(self.ceX[i]@x + self.ceB[i], self.ceY),
              (self.bmX[i]@x + self.bmB[i], self.bmY),
              ((self.lX[i]@x + self.lB[i]), self.lY)] # convert volume -> lattice constant
        r2c,r2b,r2l = [linregress(y_,y)[2]**2  for y_,y in ys]
        norm = 1. + self.bmw + self.lw
        return  (r2c + self.bmw*r2b + self.lw*r2l)/norm

    ##########
    # HELPER #
    ##########

    def metadata(self) -> dict:
        '''Summary of parameters used to determine fit inputs'''
        uid = hash_([self.calc[x] for x in ['pw','econv','decays','msb','fx']]+
                    [self.bmw,self.lw,self.lam_reg,self.cons])
        md = dict(calc   = self.calc,
                  params = dict(bmw=self.bmw,lw=self.lw,
                                reg=self.lam_reg,pre_ce=self.pre_ce,c = self.cons,
                                pre_bm=self.pre_bm,pre_l=self.pre_l),
                  uid = uid)
        return md

    @classmethod
    def data_xy(cls,
                db          : Connection,
                calc        : int,
                fx          : array,
                bmw         : float,
                lw          : float
               ) -> T[Arrs,Arrs,Arrs,Arrs,Arrs,Arrs,list,list,list,float,float,float]:
        '''
        Assemble X,Y matrices given a database connection
        Returns:
            ceXs,bmXs,lXs - lists of length 5, each having a matrix of N rows and 64 columns
            ceBs,bmBs,lBs - lists of length 5, each having a 64 element vector
            ceY,bmY,lY    - N rows 64 column matrices
            pre_c,pre_b,pre_l - lists of length 5, each having a float
        '''
        ceX_,bmX_,lX_,ceB_,bmB_,lB_,ceY,bmY,lY = [],[],[],[],[],[],[],[],[]

        for a_ce,a_bm,a_l,b_ce,b_bm,b_l,ce,bm,l,name in                         \
                [[safeLoad(x) for x in [x1,x2,x3,x4,x5,x6]]+[x7,x8,x9,x10]
                    for x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 in sqlselect(db,q3,[calc,])]:

            if a_ce and ce: ceX_.append(a_ce); ceB_.append(b_ce); ceY.append(float(ce))
            if a_bm and bm: bmX_.append(a_bm); bmB_.append(b_bm); bmY.append(float(bm))
            if a_l and l:   lX_.append(a_l);   lB_.append(b_l);   lY.append(float(l))

        Xs  = tuple(map(array,[ceX_,bmX_,lX_]))
        Bs = tuple(map(array,[ceB_,bmB_,lB_]))

        ceXs,bmXs,lXs = [[x[:,i,:].tolist()  for i in range(5)] for x in  Xs]
        ceBs,bmBs,lBs = [[x[:,i].tolist()    for i in range(5)] for x in  Bs]

        ce_avg,bm_avg,l_avg = [abs(average(x)) for x in [ceY,bmY,lY]]
        pre_c = 1 /  ce_avg
        pre_l = lw * ce_avg / l_avg
        pre_b = bmw * ce_avg / bm_avg

        return ceXs,bmXs,lXs,ceBs,bmBs,lBs,ceY,bmY,lY,pre_c,pre_b,pre_l

    def write(self, pth : str)->None:
        root = join(environ['FUNCTIONALS_ROOT'],'functionals/scripts/fit/')
        def write(fi:str,x:Any)->None:
            with open(join(pth,fi)+'.json','w') as file: dump(x,file)

        # Write to directory
        #--------------------
        write('metadata', self.metadata())
        write('data', [[x.tolist() for x in self.ceX],
                       [x.tolist() for x in  self.bmX],
                       [x.tolist() for x in self.lX],
                       [x.tolist() for x in self.ceB],
                       [x.tolist() for x in self.bmB],
                       [x.tolist() for x in self.lB],
                       self.ceY.tolist() ,self.bmY.tolist() ,self.lY.tolist() ,
                       [x.tolist() for x in self.c_A_eq],
                       [x.tolist() for x in self.c_b_eq],
                       [x.tolist() for x in self.c_A_lt],
                       [x.tolist() for x in self.c_b_lt],
                       [x.tolist() for x in self.C_A_eq],
                       [x.tolist() for x in self.C_b_eq],
                       [x.tolist() for x in self.C_A_lt],
                       [x.tolist() for x in self.C_b_lt],
                       [x.tolist() for x in self.X],
                       [y.tolist() for y in self.Y]])
        copyfile(root+'runfit.py',join(pth,'runfit.py'))
        copyfile(root+'subfit.sh',join(pth,'subfit.sh'))
        system('chmod 755 '+join(pth,'subfit.sh'))


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
