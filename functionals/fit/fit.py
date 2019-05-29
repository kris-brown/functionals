# External
from typing import (Any,
                    Set      as S,
                    Dict     as D,
                    List     as L,
                    Tuple    as T,
                    Optional as O,
                    Callable as C)

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
                            stack, array_equal, all as npall, gradient,warnings)
from numpy.linalg   import lstsq,norm,inv,LinAlgError      # type: ignore

from scipy.optimize import least_squares                   # type: ignore

from scipy.stats        import linregress         # type: ignore
from psycopg2           import connect # type: ignore
from plotly.offline     import plot # type: ignore
from plotly.graph_objs  import Scatter, Figure, Layout, Bar # type: ignore

# Internal
from functionals.fit.constraint import Constraint,all_cons
from functionals.fit.data       import Data, Datum

Connection = Any
Arrs       = L[array]
################################################################################
# CONSTANTS
###########
root = __file__.split('functionals')[0]

beef = array([
    1.18029330e+00,8.53027860e-03,-1.02312143e-01,6.85757490e-02,-6.61294786e-03,-2.84176163e-02,5.54283363e-03,3.95434277e-03,
    -1.98479086e-03,1.00339208e-01,-4.34643460e-02,-1.82177954e-02,1.62638575e-02,-8.84148272e-03,-9.57417512e-03,9.40675747e-03,
    6.37590839e-03,-8.79090772e-03,-1.50103636e-02,2.80678872e-02,-1.82911291e-02,-1.88495102e-02,1.69805915e-07,-2.76524680e-07,
    1.44642135e-03,-3.03347141e-03,2.93253041e-03,-8.45508103e-03,6.31891628e-03,-8.96771404e-03,-2.65114646e-08,5.05920757e-08,
    6.65511484e-04,1.19130546e-03,1.82906057e-03,3.39308972e-03,-7.90811707e-08,1.62238741e-07,-4.16393106e-08,5.54588743e-08,
    -1.16063796e-04,8.22139896e-04,-3.51041030e-04,8.96739466e-04,2.09603871e-08,-3.76702959e-08,2.36391411e-08,-3.38128188e-08,
    -5.54173599e-06,-5.14204676e-05,6.68980219e-09,-2.16860568e-08,9.12223751e-09,-1.38472194e-08,6.94482484e-09,-7.74224962e-09,
    7.36062570e-07,-9.40351563e-06,-2.23014657e-09,6.74910119e-09,-4.93824365e-09,8.50272392e-09,-6.91592964e-09,8.88525527e-09])

bias = array([1.3**(x%8 + x//8) for x in range(64)]) # unequally weigh matrix elements

q2 = '''SELECT ce_scale,bm_scale,lc_scale,reg,consts
        FROM fitparams WHERE fitparams_id = %s'''

q3 = '''SELECT a_ce,a_bm,a_l,b_ce,b_bm,b_l,expt_ce,expt_bm,expt_vol,name,ce
        FROM bulks JOIN job ON job=job_id
        WHERE calc = %s'''

q4 = '''SELECT pw,data,a11,a12,a13,a14,a15,msb
        FROM calc C JOIN beef B ON B.functional = C.functional
        WHERE calc_id = %s'''

################################################################################
# Helper functions
#-----------------

def sqlselect(conn : Connection, q : str, binds : list = []) -> L[tuple]:
    with conn.cursor() as cxn: # type: ignore
        cxn.execute(q,vars=binds)
        return cxn.fetchall()

def hash_(x : Any)->str: return sha512(str(x).encode()).hexdigest()[:10]

def safeLoad(x:Any)->Any: return x if x is None else loads(x)

###############################################################################
# Classes
#########

class Traj(object):
    '''
    Result of fitting: sequence from unconstrained to overconstrained solutions
    '''
    threshold = -1 # optimimum-point-finding parameter

    def __init__(self,
                 data     : Data,
                 cons     : C[[array,bool],float],
                 ces : float, bms : float, lcs : float
                 ) -> None:
        self.data  = data
        self.cons  = cons
        self.ces = ces; self.bms = bms; self.lcs = lcs

        # Initialize internal state
        self.errs = {x:[] for x in ['cost','ce','bm','lc','vol','cviol','cviol_all']} # type: D[str,L[float]]
        self.xs   = [] # type: L[array]
        self.time = time()
        # self._opt = None # type: O[int]
        self.done = False

    def __len__(self)->int:
        return len(self.xs)

    def __eq__(self,other:object)->bool:
        return False if not isinstance(other,Traj) else \
                (self.xs,self.errs) != (other.xs,other.errs)

    # Serialization/Deserialization
    def to_dict(self)->dict:
        return {**dict(x    = [x.tolist() for x in self.xs],
                       time = self.time, opt = self.opt, bms = self.bms,
                       ces = self.ces, lcs = self.lcs),**self.errs,}

    @classmethod
    def from_dict(cls,d:dict,datapth:str=None,decay:int=None)->'Traj':
        '''
        Load a completed traj instance (cannot add any more data, which would require a Fit object)
        '''
        if datapth:
            assert decay is not None
            with open(datapth,'r') as f: data = Data.from_list(load(f)[decay])
        else:
            assert decay is None
            data = Data(set())

        t = cls(data,lambda x,y:0.,d['ces'],d['bms'],d['lcs'])
        t.xs = d['x']; t.time = d['time']
        for k in ['cost','ce','bm','lc','cviol','cviol_all']:
            t.errs[k] = d[k]
        t._opt = t.getopt()
        t.done = True
        return t

    def add(self,x:array)->None:
        '''
        Add a step to the trajectory. Compute metrics upon insert
        '''
        assert not self.done
        keys = ['ce','bm','vol','lc']
        scale = dict(ce=self.ces,bm=self.bms,vol=self.lcs)
        self.xs.append(x)
        for k in keys:
            self.errs[k].append(self.data.mse(x,k))

        self.errs['cost'].append(sum([self.errs[k][-1]/scale[k] if scale[k]>0 else 0
                                      for k in keys[:3]]))
        viol    = self.cons(x,False)
        allviol = self.cons(x,True)
        self.errs['cviol'].append(viol); self.errs['cviol_all'].append(allviol)

        # COMMENT OUT LATER
        # if len(self)%500 == 1:
        #     print('\n%d'%len(self))
        #     for k in keys: print(k,self.errs[k][-1])
        #     self.resid(['lc'],step=len(self.xs)-1)
        #     import pdb;pdb.set_trace()

    def end(self)->None:
        '''Call this when an optimization is finished'''
        self.time = time() - self.time
        self._opt = self.getopt()
        self.done = True

    @property
    def x(self)->array:
        '''The "current" state of the trajectory'''
        return array(self.xs[-1])

    @property
    def opt(self)->int:
        '''
        Only compute optimum of constraint/error tradeoff once
        '''
        assert self.done; assert self._opt is not None
        return self._opt

    def getopt(self) -> int:
        '''
        Identify an optimum value (index) along the trajectory
        '''
        warnings.filterwarnings('ignore')
        try:
            cost,cviol = self.errs['cost'],self.errs['cviol']
            dic = {c : i for i,c in enumerate(cviol)} # remove duplicates
            cv, cs = self.pareto_frontier(cviol,cost)
            if len(cv)<3: return dic[cv[0]]
            dcosts = gradient(cs,cv)
            for c,dcost in zip(cv,dcosts):
                if dcost > self.threshold: break
            #import matplotlib.pyplot as plt # type: ignore
            #print(c);plt.scatter(cviol,cost); plt.show();import pdb;pdb.set_trace()
            return dic[c]
        except Exception as e:
            print(e);import pdb;pdb.set_trace();assert False

    @staticmethod
    def pareto_frontier(Xs:list, Ys:list, maxX:bool = False, maxY:bool = False) -> T[list,list]:
        '''
        By default, get lower left corner of a curve
        '''
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse = maxX)
        p_front = [myList[0]]
        for pair in myList[1:]:
            if maxY:
                if pair[1] >= p_front[-1][1]: p_front.append(pair)
            else:
                if pair[1] <= p_front[-1][1]: p_front.append(pair)
        p_frontX, p_frontY = zip(*p_front)
        return p_frontX, p_frontY

    def plot(self,met:str)->Figure:
        '''
        Plots a metric over the trajectory
        '''
        assert met in ['ce','bm','lc','vol','cost']
        xs = self.errs['cviol']
        ys = self.errs[met]
        data = [Scatter(x=xs,y=ys,name='fit',mode='markers',
                             text = [str(x) for x in range(len(xs))],
                             hoverinfo = 'text')]
        annotations =[dict(x=xs[self.opt], y=ys[self.opt],xref ='x',yref='y',text='opt')]

        layout = Layout(title= 'Trajectory of constrained optimization', hovermode= 'closest',
                        xaxis= dict(title= 'Constraint Violation', ticklen= 5, zeroline= False, gridwidth= 2),
                        yaxis= dict(title= met, ticklen= 5, gridwidth= 2,),
                        annotations = annotations)

        fig = Figure(data=data,layout=layout)
        plot(fig,filename='temp0.html')
        return Figure

    def resid(self, met:L[str]=None, step:int=None, scale : dict = None)->None:
        '''
        Plots residual BarPlot across all materials
        '''
        assert len(self.data)>0, 'need data'
        mets = met or ['ce','bm','lc']
        scales = {} if scale is None or len(mets)==1 else scale
        data = []
        s = step if step is not None else self.opt
        x = array(self.xs[s])
        for m in mets:
            mats,ys = [],[]
            for mat,d in sorted(getattr(self.data,m).items()):
                mats.append(mat); ys.append( d.err(x) / scales.get(m,1) )
            data.append(Bar(x=mats,y=ys,name=m))
        if scales:
            args = [scales.get(x,1) for x in ['ce','bm','lc']]
            scalestr = '(scaled ce/{},bm/{},lc/{})'.format(*args)
        else:
            scalestr = ''
        layout = Layout(title= 'Residuals '+scalestr, hovermode= 'closest',
                        xaxis= dict(title= 'Material'),
                        yaxis= dict(title= 'Error', ticklen= 5, gridwidth= 2,))

        fig = Figure(data=data,layout=layout)
        plot(fig,filename='temp0.html')
        return Figure


class Fit(object):
    """
    Fit a 8 x 8 matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """
    def __init__(self,
                 # DFT parameters
                 calc : dict, data : T[Data,Data,Data,Data,Data],
                 # Fit hyperparameters
                 cons : L[str], reg  : float,
                 cescale  : float, bmscale : float, lcscale : float,
                 ) -> None:
        self.data = data; self.calc = calc
        self.cons = set(cons); self.lam_reg = reg;
        self.cescale = cescale; self.bmscale = bmscale; self.lcscale = lcscale

        #print('force only CE')
        #self.cescale=0.1;self.bmscale = 0; self.lcscale = 0;

        self.c_A_eq, self.c_b_eq, self.c_A_lt, self.c_b_lt = Constraint.const_xy(calc['msb'],calc['decays'],cons)
        self.C_A_eq, self.C_b_eq, self.C_A_lt, self.C_b_lt = Constraint.const_xy(calc['msb'],calc['decays'])

        # Constant parameters for least_squares
        self.kwargs = dict(method   = 'trf',ftol = 1e-12, xtol = 1e-12, gtol = 1e-12,bounds = ([-float('inf')]*64,[float('inf')]*64),
          x_scale  = 'jac', loss     = 'linear', max_nfev = 50000, verbose  = 0) # type: dict

    def __eq__(self, other : object) -> bool:
        fields = ['calc','cons','reg','cescale','bmscale','lcscale','data']
        return False if not isinstance(other,Fit) else \
             all([getattr(self,x)==getattr(other,x) for x in fields])

    ###############
    # MAIN METHOD #
    ###############
    def cv(self, n : int = 30) -> D[str,L[float]]:
        '''
        Cross validate on 30 random 50/50 splits of the data
        '''
        keys = ['ce','bm','lc']
        results  = {k:[] for k in keys} # type: D[str,L[float]]
        for i in range(n):
            print('cv %d'%i)
            train,test = self.data[2].split()
            fitresult  = self.fit(2,train) # care about last step of i=2
            opt_x      = fitresult.xs[fitresult.opt]
            for k in keys: results[k].append(test.mse(opt_x,k))
        return results

    def fit(self, i : int, data_ : Data = None) -> Traj:
        '''
        Returns: trajectory of solutions until convergence
        '''
        data = data_ or self.data[i] # default: use all data
        t = Traj(data,lambda x,a: sum(abs(self.c_loss(x,i,a))),self.cescale,self.bmscale,self.lcscale)

        X,Y = data.xy(self.cescale,self.bmscale,self.lcscale)

        #print('adding BEEF instead')
        #t.add(beef);import pdb;pdb.set_trace()
        t.add(lstsq(X,Y,rcond=None)[0])
        dx,constr  = 0.,1e-10    # negligible constraints initially
        c_reg      = self.lam_reg/10000 # initial regularization penalty = 1/1000 of final

        # break when constraints are high AND negligible change in x
        while (constr < 1e5) or (dx > self.kwargs['xtol']):
            # Do fitting
            kwargs = dict(Xmat = X, Y = Y, constr = constr, c_reg = c_reg, i = i)
            res = least_squares(self.cost, t.x, jac = self.dcost, kwargs = kwargs, **self.kwargs)
            # Do updates
            dx      = norm(res.x - t.x)          # how much has solution changed
            c_reg  += (self.lam_reg - c_reg) / 5 # asymptotically increase regularization
            constr *= 1.01                       # increase constraints
            t.add(res.x)                         # update trajectory

        t.end()
        return t


    ################
    # CONSTRUCTORS #
    ################
    @classmethod
    def from_db(cls, db : str, fp_id : int, calc_id : int)->'Fit':
        with open(db,'r') as fi:
            kwargs = load(fi)
            kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
            conn = connect(**kwargs)

        # Extract fitting parameters
        ces,bms,lcs,reg_,constrs_ = sqlselect(conn,q2,[fp_id])[0]
        constraints = constrs_.split()
        ce_scale,bm_scale,lc_scale,reg = map(float,[ces,bms,lcs,reg_])

        # Extract calculation parameters
        pw_,fx_,a1_,a2_,a3_,a4_,a5_,msb_ = sqlselect(conn,q4,[calc_id])[0]
        pw,a11,a12,a13,a14,a15,msb = map(float,[pw_,a1_,a2_,a3_,a4_,a5_,msb_])
        calc = dict(pw=pw,decays=[a11,a12,a13,a14,a15],msb=msb,fx=loads(fx_))

        # Extract fitted data
        data = cls.db_data(db=conn,calc=calc_id)

        return cls(cons = constraints, calc=calc, cescale=ce_scale,bmscale=bm_scale,lcscale=lc_scale,
                   reg = reg,data=data)

    @classmethod
    def from_json(cls,pth : str)->'Fit':
        with open(join(pth,'metadata.json'),'r')    as fi: md = load(fi)
        with open(join(pth,'data.json'),'r')        as fi: ds = load(fi)
        p = md['params']
        d1,d2,d3,d4,d5 = [Data.from_list(d) for d in ds]
        return cls(cons=p['c'],calc=md['calc'],data=(d1,d2,d3,d4,d5),
                    reg=p['reg'],cescale=p['ce_scale'],bmscale=p['bm_scale'],lcscale=p['lc_scale'])

    #######################
    # COST/LOSS FUNCTIONS #
    #######################
    # Constraints
    @staticmethod
    def relu(x:array) -> array: return x * (x > 0)
    def c_lt_loss(self,x : array, i : int, all : bool = False) -> array:
        if all: return self.relu(self.C_A_lt[i] @ x - self.C_b_lt)
        elif len(self.c_A_lt[i]): return self.relu(self.c_A_lt[i] @ x - self.c_b_lt)
        else: return empty(0)
    def jc_lt_loss(self,x : array, i : int) -> array:
        if len(self.c_A_lt[i]): return self.c_A_lt[i] * tile(((self.c_A_lt[i] @ x - self.c_b_lt) > 0).reshape(-1,1),(1,64))
        else:                   return zeros((0,64))
    def c_eq_loss(self,x  : array, i : int, all : bool = False) -> array:
        if all: return self.C_A_eq[i] @ x - self.C_b_eq #use all known constraints
        elif len(self.c_A_eq[i]):  return self.c_A_eq[i] @ x - self.c_b_eq
        else: return empty(0)
    def jc_eq_loss(self,_ : array, i : int) -> array:
        return self.c_A_eq[i] or zeros((0,64))

    def c_loss(self,x  : array, i : int, all : bool = False) -> array:
        return concat((self.c_lt_loss(x,i,all),  self.c_eq_loss(x,i,all)))
    def dc_loss(self,x : array, i : int) -> array:
        return vstack((self.jc_lt_loss(x,i),self.jc_eq_loss(x,i)))

    # Regularization
    @staticmethod
    def reg(x  : array) -> array: return multiply(bias,x)
    @staticmethod
    def dreg(_ : array) -> array: return diag(bias)

    # Data
    def loss(self, x : array, X : array, Y : array) -> array: return X @ x - Y
    def dloss(self,X : array) -> array: return X

    # All together: data, regularization, constraints
    def cost(self, x:array,Xmat : array, Y : array, i : int,  constr : float, c_reg : float) -> array:
        return concat((self.loss(x, Xmat, Y), constr * self.c_loss(x, i), c_reg * self.reg(x)))

    def dcost(self, x:array, Xmat : array, Y : array, i : int, constr : float, c_reg : float) -> array:
        return vstack((self.dloss(Xmat), constr * self.dc_loss(x, i), c_reg * self.dreg(x)))

    ##########
    # HELPER #
    ##########

    def metadata(self) -> dict:
        '''Summary of parameters used to determine fit inputs'''
        uid = hash_([self.calc[x] for x in ['pw','decays','msb','fx']]+
                    [self.cescale,self.bmscale,self.lcscale,self.lam_reg,list(sorted(self.cons))])
        md = dict(calc   = self.calc,uid = uid,
                  params = dict(ce_scale=self.cescale,bm_scale=self.bmscale,lc_scale=self.lcscale,
                                reg=self.lam_reg,c = list(self.cons)),)
        return md

    @classmethod
    def db_data(cls,
                db    : Connection,
                calc  : int,
               ) -> T[Data,Data,Data,Data,Data]:
        '''
        Assemble datasets given a database connection
        '''
        datas = [set(),set(),set(),set(),set()] # type: L[S[Datum]]
        for a_ce,a_bm,a_lc,b_ce,b_bm,b_lc,ce,bm,lc,name,ce_calc in                         \
                [[safeLoad(x) for x in [x1,x2,x3,x4,x5,x6]]+[x7,x8,x9,x10,x11]
                    for x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11 in sqlselect(db,q3,[calc,])]:

            for i in range(5):
                if a_ce and ce and ce_calc:
                    datas[i].add(Datum(name,'ce',a_ce[i],b_ce[i],float(ce)))
                if a_bm and bm:
                    datas[i].add(Datum(name,'bm',a_bm[i],b_bm[i],float(bm)))
                if a_lc and lc:
                    datas[i].add(Datum(name,'lc',a_lc[i],b_lc[i],float(lc)))

        d1,d2,d3,d4,d5 = [Data(d) for d in datas]
        return d1,d2,d3,d4,d5

    def write(self, pth : str) -> None:
        rootpth = join(root,'functionals/functionals/scripts/fit/')
        def write(fi:str,x:Any)->None:
            with open(join(pth,fi)+'.json','w') as file: dump(x,file)

        # Write to directory
        #--------------------
        write('metadata', self.metadata())
        write('data', [d.to_list() for d in self.data])
        copyfile(rootpth+'runfit.py',join(pth,'runfit.py'))
        copyfile(rootpth+'subfit.sh',join(pth,'subfit.sh'))
        system('chmod 755 '+join(pth,'subfit.sh'))
