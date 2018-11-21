# External Modules
from typing import Callable as C, Optional as O, List as L, Tuple as T
from abc    import ABCMeta, abstractmethod
from json   import load
from numpy  import array,sum,multiply,heaviside,exp,arange # type: ignore
from plotly.graph_objs import Figure,Layout,Scatter         # type: ignore
from plotly.offline    import plot # type: ignore

# Internal
from dbgen                     import flatten
from functionals.fit.utilities import LegendreProduct
from functionals.fit           import Fit
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
        data   = flatten([p.plot() for p in ps])
        layout = Layout(title = 'Functionals',
                      xaxis = dict(title = 's'),
                      yaxis = dict(title = 'Fx'))
        fig = Figure(data = data, layout = layout)
        plot(fig,filename='temp0.html')

    def plot(self) -> L[dict]:
        ss     = arange(0.,5,0.1)
        alphas = array([0,1]) if self.mgga else [1]
        styles = ['solid','dot','dash','dashdot']

        def str2int(lst : str) -> int:
            return sum(array(list(map(ord,lst))))**3+10000000000

        color  = '#'+hex(str2int(self.name))[-6:]

        out    = []
        for sty,a in zip(styles,alphas):
            lab = self.name + r' ($\alpha$=%d)'%a if self.mgga else self.name
            out.append(Scatter(x   = ss,
                        y       = [self.apply(s,a) for s in ss],
                        mode    = 'lines',
                        name    = lab,
                        line    = dict(color   = color,
                                       dash    = sty,
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
    def __init__(self, A:array, name : O[str] = None) -> None:
        self.A = A
        self.x = A.flatten()
        self.N,self.M = A.shape
        assert self.N == self.M
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

        P = array([[LegendreProduct(s,a,j,i)
                        for j in range(self.N)]
                        for i in range(self.M)])

        return sum(multiply(self.A,P))

    def costs(self,data:str,const:str,nlconst:str)->T[float,float,float,float]:
        '''Run a single fitting step to compute costs and constraint violation'''
        f = Fit(data = data, constraints = const, nlconstraints = nlconst,
                maxit = 1, n = self.N, bound = 1, ifit = False,
                gridden = 3, bm_weight = 1, lat_weight = 1)
        classes = ['ce','bm','lat']
        r2_ce,r2_bm,r2_lat = [f.r2(self.x,c) for c in classes]
        return r2_ce,r2_bm,r2_lat,f.constr_violation(self.x)

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
    fx  = exp(-c1x*alpha/(1-alpha))*theta_1a - dx*exp(c2x/(1-alpha))*theta_a1
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
with open('/Users/ksb/functionals/data/beef.json') as f:
    beefcoeff = array(load(f))

BEEF = FromMatrix(beefcoeff.reshape(8,8),'beef')
