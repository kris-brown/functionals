# External
from typing import (Any,
                    Set      as S,
                    Dict     as D,
                    List     as L,
                    Union    as U,
                    Tuple    as T,
                    Optional as O,
                    Callable as C)

from os.path      import join, exists
from shutil       import copyfile
from os           import system, listdir, getcwd
from json         import dumps, loads, dump, load
from hashlib      import sha512
from collections  import defaultdict
from tempfile     import TemporaryDirectory
from numpy          import (empty,tile,sum,ones,array,linspace,  # type: ignore
                            zeros,inf,vstack,concatenate as concat,diag, multiply,
                            mean,minimum,maximum,exp, sign, diff, where)
from numpy.linalg   import lstsq,norm      # type: ignore

from sklearn            import neighbors                    # type: ignore
from scipy.optimize     import least_squares                # type: ignore
from scipy.interpolate  import UnivariateSpline             # type: ignore
from psycopg2           import connect                      # type: ignore
from plotly.offline     import plot                         # type: ignore
from plotly.graph_objs  import Scatter, Figure, Layout, Bar # type: ignore

# Internal
from functionals.fit.constraint import Constraint,ConstraintType,all_cons
from functionals.fit.data       import Data, Datum
from functionals.fit.functional import FromMatrix,Functional,PBE,BEEF,SCAN
Connection = Any
Arrs       = L[array]
Metric     = U[str,D[str,float]]
################################################################################
# CONSTANTS
###########
errtypes = ['ce','bm','lc']

root = '/'+join(*__file__.split('/')[:-3])

beef = array([1.18029330e+00,8.53027860e-03,-1.02312143e-01,6.85757490e-02,-6.61294786e-03,-2.84176163e-02,5.54283363e-03,3.95434277e-03,-1.98479086e-03,1.00339208e-01,-4.34643460e-02,-1.82177954e-02,1.62638575e-02,-8.84148272e-03,-9.57417512e-03,9.40675747e-03,6.37590839e-03,-8.79090772e-03,-1.50103636e-02,2.80678872e-02,-1.82911291e-02,-1.88495102e-02,1.69805915e-07,-2.76524680e-07,1.44642135e-03,-3.03347141e-03,2.93253041e-03,-8.45508103e-03,6.31891628e-03,-8.96771404e-03,-2.65114646e-08,5.05920757e-08,6.65511484e-04,1.19130546e-03,1.82906057e-03,3.39308972e-03,-7.90811707e-08,1.62238741e-07,-4.16393106e-08,5.54588743e-08,-1.16063796e-04,8.22139896e-04,-3.51041030e-04,8.96739466e-04,2.09603871e-08,-3.76702959e-08,2.36391411e-08,-3.38128188e-08,-5.54173599e-06,-5.14204676e-05,6.68980219e-09,-2.16860568e-08,9.12223751e-09,-1.38472194e-08,6.94482484e-09,-7.74224962e-09,7.36062570e-07,-9.40351563e-06,-2.23014657e-09,6.74910119e-09,-4.93824365e-09,8.50272392e-09,-6.91592964e-09,8.88525527e-09])

bias = array([1.2**(x%8 + x//8) for x in range(64)]) # unequally weigh matrix elements

q3 = '''SELECT a_ce,a_bm,a_l,b_ce,b_bm,b_l,expt_ce,expt_bm,expt_vol,name,ce
        FROM bulks JOIN job ON job=job_id
        WHERE calc = %s'''

q4 = '''SELECT pw,data,fitdata
        FROM calc JOIN functional on functional=functional_id
        WHERE calc_id = %s'''

################################################################################
# Helper functions
#-----------------
def pltstr(x:Figure,show:bool=True)->str:
    with TemporaryDirectory() as tmpdirname:
        plot(x,filename=tmpdirname+'/temp0.html',auto_open=show)
        with open(tmpdirname+'/temp0.html','r') as f: out = f.read()
    return out

def sqlselect(conn : Connection, q : str, binds : list = []) -> L[tuple]:
    with conn.cursor() as cxn: # type: ignore
        cxn.execute(q,vars=binds)
        return cxn.fetchall()

def hash_(x : Any)->str: return sha512(str(x).encode()).hexdigest()[:10]

def safeLoad(x:Any)->Any: return x if x is None else loads(x)

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

###############################################################################
# Classes
#########

class Traj(object):
    '''
    Result of fitting: sequence from unconstrained to overconstrained solutions.
    '''
    BOUNDARY_THRESHOLD = 0.9
    BEND_THRESHOLD = 5

    def __init__(self, train : Data, test : Data,
                 cons : S[Constraint]) -> None:
        self.train = train
        self.test  = test
        self.cons  = cons
        # Initialize state
        self.xs        = [] # type: L[array]
        self.train_err = defaultdict(list) # type: D[str,L[float]]
        self.test_err  = defaultdict(list) # type: D[str,L[float]]
        self.cviol     = defaultdict(list) # type: D[str,L[float]]
        self.frac      = [] # type: L[float]
        self.bends     = [] # type: L[int]

    def __len__(self)->int: return len(self.xs)

    def __eq__(self,other:object)->bool:
        return False if not isinstance(other,Traj) else \
                (vars(self)) != (vars(other))

    # Serialization/Deserialization
    def to_dict(self)->dict:
        return dict(train = self.train.to_list(),test = self.test.to_list(),
                    cons  = [c.to_dict() for c in self.cons],
                    xs    = [x.tolist() for x in self.xs],
                    **{x:getattr(self,x) for x in ['train_err','test_err','cviol','frac','bends']})

    @classmethod
    def from_dict(cls, d : dict, dx : int = 1) -> 'Traj':
        '''Inverse of to_dict.'''
        t = cls(train = Data.from_list(d['train']),test = Data.from_list(d['test']),
                cons ={Constraint.from_dict(x) for x in d['cons']})
        t.xs = [array(x) for x in d['xs'][::dx]]
        t.frac, t.bends = d['frac'], d['bends']

        for x in ['train_err','test_err','cviol']: setattr(t,x,{k:v[::dx] for k,v in d[x].items()})

        return t

    def add(self, x : array) -> None:
        '''Add a step to the trajectory. Compute metrics upon insert.'''
        self.xs.append(x)
        for k in errtypes:
            self.train_err[k].append(self.train.mse(x,k))
            self.test_err[k].append(self.test.mse(x,k))

        for c in self.cons:
            self.cviol[c.name].append(c.err(x))

        fx = FromMatrix(x,a1=4.9479,msb=1)
        vals = [[fx.apply(s,a) for s in linspace(0,5,20)] for a in [0,1,]]
        fracs = [float(mean([v > 0 and v < 1.804 for v in val])) for val in vals]
        self.frac.append(sum(fracs)/2)
        dirs = [sign(diff(val)) for val in vals]
        bends= [len(where(diff(d))[0]) for d in dirs]
        self.bends.append(sum(bends))

    def gap(self, met : Metric)->array:
        '''Gap between train and test'''
        return maximum(self.weighted_err(met,train = False) - self.weighted_err(met,train = True),0.0)

    def weighted_err(self,met:Metric,train:bool) -> array:
        metdic = {met:1.} if isinstance(met,str) else met
        it = self.train_err if train else self.test_err
        return sum([metdic.get(k,0)*array(v) for k,v in it.items()],axis=0)

    def cviols(self,cons:L[str]=None)->array:
        '''Sum all constraint violations found in cons'''
        return sum([array(x) for k,x in self.cviol.items() if (cons is None) or (k in cons)],axis=0)

    def opt_err(self) -> T[float,float,float]:
        a,b,c = [self.test_err[x][self.opt] for x in errtypes]
        return a,b,c

    @property
    def opt(self) -> int:
        for i, (f, b) in enumerate(zip(self.frac,self.bends)):
            if f > self.BOUNDARY_THRESHOLD and b < self.BEND_THRESHOLD:
                return i
        import pdb;pdb.set_trace()
        raise ValueError()
    #################
    # VISUALIZATION #
    #################
    def fxviz(self,name:str='fit',show:bool=True, opt:int=None) -> str:
        fx = FromMatrix(A=array(self.xs[opt or self.opt]), a1=4.9479, msb=1., name=name)
        fig = Functional.plots([fx,PBE,BEEF,SCAN])
        return pltstr(fig,show=show)

    def plot(self,met:Metric,cons:L[str]=None,show:bool=True) -> str:
        '''Plots a metric over the trajectory.'''
        xs = self.cviols(cons)

        trainys = self.weighted_err(met,train=True)
        data = [Scatter(x=xs,y=trainys,name='train',mode='markers',
                             text = [str(x) for x in range(len(xs))],
                             hoverinfo = 'text')]

        if self.test != self.train:
            testys  = self.weighted_err(met,train=False)
            data.append(Scatter(x=xs,y=testys,name='test',mode='markers',
                                     text = [str(x) for x in range(len(xs))],
                                     hoverinfo = 'text'))
        annotations =[]  # type: list

        layout = Layout(title= 'Trajectory of constrained optimization', hovermode= 'closest',
                        xaxis= dict(title= 'Constraint Violation', ticklen= 5, zeroline= False, gridwidth= 2),
                        yaxis= dict(title= str(met), ticklen= 5, gridwidth= 2,),
                        annotations = annotations)

        fig = Figure(data=data,layout=layout)

        return pltstr(fig,show=show)

    def resid(self, step:int=-1, train : bool = False, show:bool=True) -> str:
        '''Plots residual BarPlot across all materials.'''
        data     = self.train if train else self.test
        datalist = []
        x        = array(self.xs[step])
        for m in errtypes:
            mats,ys = [],[]
            for d in getattr(data,m):
                mats.append(d.mat); ys.append( d.err(x) )
            datalist.append(Bar(x=mats,y=ys,name=m))

        layout = Layout(title= 'Residuals ', hovermode= 'closest',
                        xaxis= dict(title= 'Material'),
                        yaxis= dict(title= 'Error', ticklen= 5, gridwidth= 2,))
        fig = Figure(data=datalist,layout=layout)
        return pltstr(fig,show=show)

class FitResult(object):
    def __init__(self,full:Traj,cv:L[Traj])->None:
        self.full = full
        self.cv   = cv

    @classmethod
    def from_pth(cls, pth : str, dx : int = 1) -> 'FitResult':
        cv = []
        try:
            with open(join(pth,'fit.json'),'r') as f:
                data = Traj.from_dict(load(f),dx=dx)
        except Exception as e:
            print(pth,e)
            import pdb;pdb.set_trace()
        for x in filter(lambda x:x[0]=='x',listdir(pth)):
            with open(join(pth,x),'r') as f:
                try:
                    cv.append(Traj.from_dict(load(f),dx=dx))
                except:
                    print('error loading x ',pth,x)
        return cls(data,cv)

    def plot_cv(self,met : Metric,show:bool=True) -> str:
        '''
        Plots a metric over the trajectory
        '''
        kwargs = dict(mode='markers',hoverinfo = 'text')
        data = []#Scatter(name='full',x=fullx,y=fully,**kwargs)]
        for i,t in enumerate(self.cv):
            x = t.cviols()
            y = t.weighted_err(met,train=True)
            z = t.weighted_err(met,train=False)
            text= [str(x) for x in range(len(x))]
            data.extend([Scatter(x=x,y=y,name='train_%d'%i,text = text,**kwargs),
                         Scatter(x=x,y=z,name='test_%d'%i,**kwargs),])

        layout = Layout(title= 'Trajectory of constrained optimization', hovermode= 'closest',
                        xaxis= dict(title= '∆Train/Test', ticklen= 5, zeroline= False, gridwidth= 2),
                        yaxis= dict(title= 'Test error', ticklen= 5, gridwidth= 2,),)

        fig = Figure(data=data,layout=layout)
        return pltstr(fig,show=show)

    def plot_transfer(self, met : Metric,cons:L[str] = None, show : bool = True) -> str:
        fullx = self.full.cviols()
        fully = self.full.weighted_err(met,train=False)
        tx,ty,tx_,ty_ = self.transferability(met,cons)

        kwargs = dict(mode='markers',text = [str(x) for x in range(len(fullx))],hoverinfo = 'text')
        data = [Scatter(name='full',x=fullx,y=fully,**kwargs),
                Scatter(name='∆Err',x=tx,y=ty,yaxis='y2',**kwargs),
                Scatter(name='∆Err_smooth',x=tx_,y=ty_,yaxis='y2',**kwargs),]
        layout = Layout(title= 'Transferability', hovermode= 'closest',
                        xaxis= dict(title= 'Constraint Violation', ticklen= 5, zeroline= False, gridwidth= 2),
                        yaxis= dict(title= str(met), ticklen= 5, gridwidth= 2,),
                        yaxis2=dict(title='∆Err',overlaying='y',side='right'))

        fig = Figure(data=data,layout=layout)
        return pltstr(fig,show=show)

    def transferability(self, met : Metric, cons : L[str]=None) -> T[array,array,array,array]:
        '''For each data point, compute average ErrTest / ErrTrain ratio'''
        xs,tsterr,trnerr = [],[],[] # type: ignore
        for traj in self.cv:
            xs.extend(traj.cviols(cons));
            tsterr.extend(traj.weighted_err(met,train = False));
            trnerr.extend(traj.weighted_err(met,train = True))

        # KNN regression
        knn = neighbors.KNeighborsRegressor(20)
        ys = maximum(array(tsterr)-array(trnerr),0.0) # minimum(array(trnerr)/array(tsterr),1.0)
        x_ = array(xs).reshape((-1,1))

        outx = self.full.cviols() #linspace(min(xs),max(xs),200)

        outy = knn.fit(X=x_, y=ys).predict(outx.reshape((-1,1)))

        #sx,sy = zip(*sorted(zip(x_,ys_))) # sort for spline input
        #f = UnivariateSpline(x = sx,y =sy, k=1)
        #outy = array([f(x) for x in outx])
        return xs,ys,outx,outy


class Fit(object):
    """
    Fit a 8 x 8 matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """
    def __init__(self,
                 # DFT parameters
                 calc : dict, data : Data,
                 # Fit hyperparameters
                 cons : L[str], reg  : float,
                 cescale  : float, bmscale : float, lcscale : float,
                 mag_scale : float,
                 ) -> None:
        self.data = data; self.calc = calc
        self.cons = set(cons); self.lam_reg = reg;
        self.cescale = cescale; self.bmscale = bmscale; self.lcscale = lcscale
        self.mag_scale=mag_scale

        self.c_A_eq, self.c_b_eq, self.c_A_lt, self.c_b_lt = ConstraintType.const_xy(cons)
        self.C_A_eq, self.C_b_eq, self.C_A_lt, self.C_b_lt = ConstraintType.const_xy()

        # Constant parameters for least_squares
        self.kwargs = dict(method = 'trf',ftol = 1e-12, xtol = 1e-12, gtol = 1e-12,
                            bounds = ([-float('inf')]*64,[float('inf')]*64),
                            x_scale = 'jac', loss = 'linear', max_nfev = 50000,
                            verbose  = 0) # type: dict

    def __eq__(self, other : object) -> bool:
        fields = ['calc','cons','reg','cescale','bmscale','lcscale','data','mag_scale']
        return False if not isinstance(other,Fit) else \
             all([getattr(self,x)==getattr(other,x) for x in fields])

    ###############
    # MAIN METHOD #
    ###############
    def fit(self, pth : str, n : int = 30) -> None:
        ''' Cross validate on 30 random 50/50 splits of the data.'''
        tpth = join(pth,'fit.json')
        if not exists(tpth):
            t = self.fittraj().to_dict()
            with open(tpth,'w') as f: dump(t,f)

        # Do cross validation
        for i in range(10):
            tpth = join(pth,'x%d.json'%i)
            trn,tst = self.data.split(i)
            if not exists(tpth):
                t = self.fittraj(train=trn, test=tst).to_dict()
                with open(tpth,'w') as f: dump(t,f)


    def fittraj(self, train : Data = None, test : Data = None) -> Traj:
        '''
        Returns: trajectory of constrained solutions to Ax=b until convergence
        '''
        if train is None and test is None:
            train,test = (self.data,self.data)
        assert train is not None; assert test is not None

        # initialize empty trajectory
        #############################
        cons  = {all_cons[k].mkcon() for k in self.cons}
        t = Traj(train, test, cons = cons,)

        # Initialize fitting
        #####################
            #print('adding BEEF instead - MAKE SURE RESIDUALS MATCH DATABASE');x=beef;t.add(x);import pdb;pdb.set_trace()

        X,Y        = train.xy(self.cescale, self.bmscale,
                              self.lcscale, self.mag_scale) # created weighted Ax=b
        x          = lstsq(X,Y,rcond=None)[0]                           # unconstrained, unregularized solution
        constr     = 1e-10                                              # negligible constraints initially
        c_reg      = self.lam_reg/10000                                 # negligible regularization initially
        dx,counter = 0., 0

        # break when constraints are high AND negligible change in x
        ############################################################
        while (constr < 1e5) or (dx > 1e-9):
            # Do fitting
            kwargs = dict(Xmat = X, Y = Y, constr = constr, c_reg = c_reg)
            res = least_squares(self.cost, x, jac = self.dcost, kwargs = kwargs, **self.kwargs)
            # Do updates
            dreg    = self.lam_reg - c_reg    # difference between regularization term and final value
            if dreg > 1e-5:                   # We are in "Startup" - regularization is changing, not constraints
                c_reg  += (dreg) / 5          # asymptotically increase regularization
            else:
                dx       = norm(res.x - x)    # how much has solution changed
                constr  *= 1.02               # increase constraints
                x        = res.x              # update current vector
                counter +=1
                if counter % 5 == 0: t.add(x) # update trajectory

        return t

    ################
    # CONSTRUCTORS #
    ################
    @classmethod
    def from_db(cls, db : str, calc_id : int, constr:str, reg:float,
                ce:float,bm:float,lc:float,mc:float) -> 'Fit':
        with open(db,'r') as fi:
            kwargs = load(fi)
            kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
            conn = connect(**kwargs)

        # Extract fitting parameters
        constraints = constr.split()
        ce_scale,bm_scale,lc_scale,reg,mag_scale = map(float,[ce,bm,lc,reg,mc])

        # Extract calculation parameters
        pw_,fx_,data_ = sqlselect(conn,q4,[calc_id])[0]
        pw = float(pw_)
        calc = dict(pw=pw,fx=loads(fx_))

        # Extract fitted data
        data = Data.from_list(loads(data_),full=True)

        return cls(cons = constraints, calc=calc, cescale=ce_scale, bmscale=bm_scale,
                   lcscale=lc_scale, mag_scale=mag_scale,reg = reg, data=data)

    @classmethod
    def from_json(cls,pth : str)->'Fit':
        with open(join(pth,'metadata.json'),'r') as fi: md = load(fi)
        with open(join(pth,'data.json'),'r')     as fi: data = load(fi)
        p = md['params']
        return cls(cons=p['c'],calc=md['calc'],data=Data.from_list(data),
                    reg=p['reg'],cescale=p['ce_scale'],bmscale=p['bm_scale'],
                    lcscale=p['lc_scale'],mag_scale=p['mag_scale'])

    #######################
    # COST/LOSS FUNCTIONS #
    #######################
    # Constraints
    @staticmethod
    def relu(x:array) -> array: return x * (x > 0)
    def c_lt_loss(self,x : array, all : bool = False) -> array:
        if all: return self.relu(self.C_A_lt @ x - self.C_b_lt)
        elif len(self.c_A_lt): return self.relu(self.c_A_lt @ x - self.c_b_lt)
        else: return empty(0)
    def jc_lt_loss(self,x : array) -> array:
        if len(self.c_A_lt): return self.c_A_lt * tile(((self.c_A_lt @ x - self.c_b_lt) > 0).reshape(-1,1),(1,64))
        else:                   return zeros((0,64))
    def c_eq_loss(self,x  : array, all : bool = False) -> array:
        if all: return self.C_A_eq @ x - self.C_b_eq #use all known constraints
        elif len(self.c_A_eq):  return self.c_A_eq @ x - self.c_b_eq
        else: return empty(0)
    def jc_eq_loss(self,_ : array) -> array:
        return self.c_A_eq or zeros((0,64))

    def c_loss(self,x  : array, all : bool = False) -> array:
        return concat((self.c_lt_loss(x,all),  self.c_eq_loss(x,all)))
    def dc_loss(self,x : array) -> array:
        return vstack((self.jc_lt_loss(x),self.jc_eq_loss(x)))

    # Regularization
    @staticmethod
    def reg(x  : array) -> array: return multiply(bias,x)
    @staticmethod
    def dreg(_ : array) -> array: return diag(bias)

    # Data
    def loss(self, x : array, X : array, Y : array) -> array: return X @ x - Y
    def dloss(self,X : array) -> array: return X

    # All together: data, regularization, constraints
    def cost(self, x:array,Xmat : array, Y : array,  constr : float, c_reg : float) -> array:
        return concat((self.loss(x, Xmat, Y), constr * self.c_loss(x), c_reg * self.reg(x)))

    def dcost(self, x:array, Xmat : array, Y : array, constr : float, c_reg : float) -> array:
        return vstack((self.dloss(Xmat), constr * self.dc_loss(x), c_reg * self.dreg(x)))

    ##########
    # HELPER #
    ##########

    def metadata(self) -> dict:
        '''Summary of parameters used to determine fit inputs'''
        uid = hash_([self.calc[x] for x in ['pw','fx']]+
                    [self.cescale,self.bmscale,self.lcscale,self.mag_scale,
                    self.lam_reg,list(sorted(self.cons))])
        md = dict(calc   = self.calc,uid = uid,
                  params = dict(ce_scale=self.cescale,bm_scale=self.bmscale,
                                lc_scale=self.lcscale, mag_scale=self.mag_scale,
                                reg=self.lam_reg,c = list(self.cons)),)
        return md

    def write(self, pth : str) -> None:
        rootpth = join(root,'functionals/scripts/fit/')
        def write(fi:str,x:Any)->None:
            with open(join(pth,fi)+'.json','w') as file: dump(x,file)

        # Write to directory
        #--------------------
        write('metadata', self.metadata())
        write('data', self.data.to_list())
        copyfile(rootpth+'runfit.py',join(pth,'runfit.py'))
        copyfile(rootpth+'subfit.sh',join(pth,'subfit.sh'))
        system('chmod 755 '+join(pth,'subfit.sh'))

# if __name__=='__main__':
#     f = FitResult.from_pth(getcwd(),dx=1)
#     # f.plot_cv('bm')
#     # f.plot_transfer('bm')
#     f.cv[0].resid(train=True)
