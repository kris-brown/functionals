# External Modules
from typing import Callable as C, List as L, Dict as D, Tuple as T
from abc    import ABCMeta, abstractmethod
from json   import load, loads
from os     import environ
from hashlib import sha224
from os.path import join
from numpy  import inf,vstack,array,sum,multiply,heaviside,exp,arange,concatenate as concat # type: ignore
from plotly.graph_objs import Figure,Layout,Scatter        # type: ignore
from plotly.offline    import plot # type: ignore
from sklearn.linear_model import LinearRegression as LinReg          # type: ignore

# Internal
from functionals.fit.utilities import LegendreProduct
from functionals.fit.fit           import Fit,sqlselect
binary = C[[float,float],float]
################################################################################
################################################################################
################################################################################
# Functional constants
#---------------------
mu_pbe   = 0.2195149727645171
kappa    = 0.804
mu       = 10./81
################################################################################
################################################################################
################################################################################
def flatten(lol:list)->list: return [item for sublist in lol for item in sublist]


# Functional classes
#-------------------
class Functional(object, metaclass = ABCMeta):
    """
    Something that implements a function (s,⍺)->Fx :: (ℝ+,ℝ+)->ℝ+
    """
    @property
    @abstractmethod
    def name(self) -> str: pass

    @property
    @abstractmethod
    def mgga(self) -> bool: pass

    @abstractmethod
    def apply(self, s : float, a : float) -> float:
        raise NotImplementedError

    @staticmethod
    def plots(ps : L['Functional']) -> None:
        assert len(ps)<7
        cs     = colors = ['rgb(255,0,0)', 'rgb(0,255,0)', 'rgb(0,0,255)', 'rgb(153,0,153)', 'rgb(0,0,0)', 'rgb(255,255,0)']
        data   = flatten([p.plot(color=c) for p,c in zip(ps,cs)])
        layout = Layout(title = 'Functionals',
                      xaxis = dict(title = 's'),
                      yaxis = dict(title = 'Fx'))
        fig = Figure(data = data, layout = layout)
        plot(fig,filename='temp0.html')

    def plot(self,color:str) -> L[dict]:
        ss     = arange(0.,5,0.1)
        alphas = array([0,1]) if self.mgga else [1]
        styles = ['solid','dot','dash','dashdot']


        out    = []
        for sty,a in zip(styles,alphas):
            lab = self.name + r' ($\alpha$=%d)'%a if self.mgga else self.name
            out.append(Scatter(x   = ss,
                        y       = [self.apply(s,a) for s in ss],
                        mode    = 'lines',
                        name    = lab,
                        line    = dict(dash    = sty,
                                       color   = color,
                                       shape   = 'spline')))
        return out

    @staticmethod
    def diff(f1:'Functional',f2:'Functional',s:L[float],a:L[float])->float:
        '''
        Return the average distance between Fx for a grid of points
        '''
        diffs = []
        for s_ in s:
            for a_ in a:
                diffs.append(abs(f1.apply(s_,a_)-f2.apply(s_,a_)))
        return sum(diffs)/len(diffs)

# Different ways of creating a functional
class FromFunc(Functional):
    """ Directly specify the enhancement factor formula """
    def __init__(self, name : str, f : binary, mgga : bool = True) -> None:
        self._name = name; self.f = f; self._mgga = mgga

    @property
    def name(self) -> str: return self._name

    @property
    def mgga(self) -> bool: return self._mgga

    def apply(self, s : float, a : float) -> float:
        return self.f(s,a)

class FromMatrix(Functional):
    """
    Implicitly define Fx formula in terms of a 2D matrix corresponding to
    coefficients for Legendre Polynomials
    """
    def __init__(self, A:array, a1 : float, msb : float, name : str = None) -> None:
        self.A = A.reshape((8,8))
        self.x = A.flatten()
        self.a1 = a1
        self.msb = msb
        assert self.A.shape == (8,8)
        self._name = name or '<no name>'

    @property
    def name(self) -> str: return self._name

    @property
    def mgga(self) -> bool: return True

    def apply(self, s : float, a : float) -> float:
        """
        A - an MxN matrix with rows corresponding to s basis functions coefficients
            and columns corresponding to alpha basis function coefficients
        Eq # 5 in mBEEF paper
        """
        P = array([LegendreProduct(s,a,self.a1,self.msb,i,j)
                                   for i in range(8) for j in range(8)])
        return sum(self.x @ P)

    @staticmethod
    def mkPlots(a1:float,msb:float,**ps : str) -> None:
        Functional.plots([FromMatrix(loads(p) if isinstance(p,str) else p,a1,msb,n)
                            for n,p in ps.items()])

    # def costs(self,dbpth:str,calc:int,decay:int) -> T[float,float,float,float]:
    #     '''Returns r2_ce, r2_bm, r2_lat, and c_viol'''
    #     with open(dbpth,'r') as f: conn = connect(**load(f),  autocommit  = True)
    #
    #     constraints = {n:1 for (n,) in sqlselect(conn,'SELECT name FROM const')}
    #     fit = Fit.from_db(db=dbpth,constraints=constraints,calc_id=calc,decay=decay,dataconstr='.',bmw=0,lw=0,reg=0,cd=0)
    #     ceX,bmX,lX,ceY,bmY,lY = fit.data
    #     p = fit.metadata()['params']
    #     c_A_eq,c_b_eq,c_A_lt,c_b_lt = fit.const_xy()
    #
    #     def relu(x : array) -> array: return x * (x > 0)
    #     def c_lt_loss(x : array) -> float: return relu(c_A_lt @ x - c_b_lt)
    #     def c_eq_loss(x  : array) -> float: return c_A_eq @ x - c_b_eq
    #     def c_loss(x : array) -> float: return concat((c_lt_loss(x), c_eq_loss(x)))
    #
    #     rc,rb,rl =  [LinReg().fit(a @ self.x, b).score(a @ self.x, b) for a,b in
    #                     [(ceX,ceY),(bmX,bmY),(lX,lY)]]
    #     return rc,rb,rl,c_loss(self.x)

    # def spaghetti(self)->np.array:
    #     #Effective number of degrees of freedom gives a scale of the ensemble
    #     #temperature T something to play with, controls the Spaghetti spread
    #     #----------------------------------------------------------------------
    #     Ndeg = 2.5                          # Equation 8 (for actual formula)
    #     min_cost = 0.5 * np.sum(self.functional().resid(self.data)[0]**2)     # Equation 9
    #     T   = 2. * min_cost / Ndeg          # Equation 11b
    #     A   = np.dot(self.data.x.T,self.data.x) / T
    #     l,U = np.linalg.eigh(A)
    #     L   = np.diag(1./np.sqrt(l))
    #     M   = np.dot(U,L)  # ensemble information
    #     return M

################################################################################
################################################################################
################################################################################
# Functional instances
#---------------------

# Functionals generated by Explicit enhancement factor functions
#----------------------------------------------------------------
PBE    = FromFunc('PBE',   lambda s,_: 1. + kappa*(1. - 1./(1. + mu_pbe*s**2 / kappa)),mgga=False)
RPBE   = FromFunc('RPBE',  lambda s,_: 1. + kappa*(1. - exp(-mu_pbe*s**2 / kappa)),mgga=False)
PBEsol = FromFunc('PBEsol',lambda s,_: 1. + kappa*(1. - 1./(1.+mu*s**2 / kappa)),mgga=False)

def fxSCAN(s:float,alpha:float)->float:
    # Scan-specific constants
    h0x,c1x,c2x,dx,b3,k1,a1 = 1.174,0.667,0.8,1.24,0.5,0.065,4.9479
    b2 = (5913/405000)**0.5
    b1 = (511/13500) / (2*b2)
    b4 = mu**2/k1 - 1606/18225 - b1**2
    # Edge conditions with numerical instability
    assert s >= 0 and alpha >= 0
    if s < 0.01: s = 0.01
    if abs(alpha-1) < 0.01: alpha = 1.001
    # Intermediate values
    theta_1a = float(heaviside(1-alpha,0.5))
    theta_a1 = float(heaviside(alpha-1,0.5))
    s2  = s**2
    x   = mu*s2*(1 + b4*s2/mu) + (b1*s2 + b2*(1-alpha)*exp(-b3*(1-alpha)**2))**2
    h1x = 1 + k1 - k1/(1+x/k1)
    gx  = 1 - exp(-a1*s**(-0.5))
    fx  = exp(-c1x * alpha/(1-alpha)) * theta_1a - dx*exp(c2x/(1-alpha))*theta_a1
    # Main output
    return (h1x + fx*(h0x-h1x))*gx

def fxSCAN_(s:float,alpha:float)->float:
    '''Comment out taylor expansion terms in gx as desired'''
    # Scan-specific constants
    h0x,c1x,c2x,dx,b3,k1,a1 = 1.174,0.667,0.8,1.24,0.5,0.065,4.9479
    b2 = (5913/405000)**0.5
    b1 = (511/13500) / (2*b2)
    b4 = mu**2/k1 - 1606/18225 - b1**2
    # Edge conditions with numerical instability
    assert s >= 0 and alpha >= 0
    if s < 0.01: s = 0.01
    if abs(alpha-1) < 0.01: alpha = 1.001
    # Intermediate values
    theta_1a = float(heaviside(1-alpha,0.5))
    theta_a1 = float(heaviside(alpha-1,0.5))
    s2  = s**2
    x   = mu*s2*(1 + b4*s2/mu) + (b1*s2 + b2*(1-alpha)*exp(-b3*(1-alpha)**2))**2
    h1x = 1 + k1 - k1/(1+x/k1)
    #gx  = 1 - exp(-a1*s**(-0.5))
    ea = exp(-a1)
    gx   = ((1-ea) - 0.5*(a1*ea)*(s-1) + a1*ea*(s-1)**2 *(0.375-0.125*a1)
            +a1 * (-0.0208333 * a1**2 + 0.1875 * a1 - 0.3125) * ea * (s - 1)**3)
    fx  = exp(-c1x * alpha/(1-alpha)) * theta_1a - dx*exp(c2x/(1-alpha))*theta_a1
    # Main output
    return (h1x + fx*(h0x-h1x))*gx

SCAN = FromFunc('SCAN',fxSCAN)

def fxMS2(s:float,alpha:float)->float:
    k,c,b = 0.504,0.14601,4.0
    p     = s**2
    F1x   = 1 + k - k/(1+mu*p/k)
    F0x   = 1 + k - k/(1+(mu*p+c)/k)
    f     = (1-alpha**2) / (1 + alpha**3 + b*alpha**6)
    return F1x + f*(F0x-F1x)

MS2 = FromFunc('MS2',fxMS2)

# From data
#----------
with open(join(environ['FUNCTIONALS_ROOT'],'data/beef.json'),'r') as f:
    beefcoeff = array(load(f))

BEEF = FromMatrix(beefcoeff.reshape(8,8),a1=float('inf'),msb=1,name='beef')
