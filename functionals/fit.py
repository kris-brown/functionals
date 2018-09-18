# External
from typing     import Tuple,List,Optional,Callable,Dict
from sys        import argv,exit
from csv        import DictReader,reader,writer
from pprint     import pformat
from operator   import add
from collections import defaultdict
from ast        import literal_eval
from math       import exp
from os.path    import exists
from numpy      import array,heaviside  # type: ignore
from matplotlib import pyplot as plt   # type: ignore
from matplotlib.axes import Axes               # type: ignore
import numpy as np                      # type: ignore
import pdb
###############################################################################
from matplotlib import rc
#rc('text', usetex=True)
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
plt.rcParams.update(pgf_with_rc_fonts)

np.set_printoptions(precision=2,linewidth=120,floatmode='fixed',suppress=True)

############
# Constants
###########

# Type Synonyms
csv_out_type = Tuple[List[str],array,array,array,array,array,array]
binary       = Callable[[float,float],float]

# BEEF
beefcoeff = array([
    1.18029330e+00,8.53027860e-03,-1.02312143e-01,6.85757490e-02,-6.61294786e-03,-2.84176163e-02,5.54283363e-03,3.95434277e-03,
    -1.98479086e-03,1.00339208e-01,-4.34643460e-02,-1.82177954e-02,1.62638575e-02,-8.84148272e-03,-9.57417512e-03,9.40675747e-03,
    6.37590839e-03,-8.79090772e-03,-1.50103636e-02,2.80678872e-02,-1.82911291e-02,-1.88495102e-02,1.69805915e-07,-2.76524680e-07,
    1.44642135e-03,-3.03347141e-03,2.93253041e-03,-8.45508103e-03,6.31891628e-03,-8.96771404e-03,-2.65114646e-08,5.05920757e-08,
    6.65511484e-04,1.19130546e-03,1.82906057e-03,3.39308972e-03,-7.90811707e-08,1.62238741e-07,-4.16393106e-08,5.54588743e-08,
    -1.16063796e-04,8.22139896e-04,-3.51041030e-04,8.96739466e-04,2.09603871e-08,-3.76702959e-08,2.36391411e-08,-3.38128188e-08,
    -5.54173599e-06,-5.14204676e-05,6.68980219e-09,-2.16860568e-08,9.12223751e-09,-1.38472194e-08,6.94482484e-09,-7.74224962e-09,
    7.36062570e-07,-9.40351563e-06,-2.23014657e-09,6.74910119e-09,-4.93824365e-09,8.50272392e-09,-6.91592964e-09,8.88525527e-09])

# Functional constants
mu_pbe   = 0.2195149727645171
kappa    = 0.804
mu       = 10./81


##############
# Functionals
#-------------
def fxPBE(s:float,_:float)->float:
    return 1. + kappa*(1. - 1./(1. + mu_pbe*s**2 / kappa))

def fxRPBE(s:float,_:float)->float:
    return 1. + kappa*(1. - np.exp(-mu_pbe*s**2 / kappa))

def fxPBEsol(s:float,_:float)->float:
    return 1. + kappa*(1.- 1./(1.+mu*s**2/kappa))

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

def fxMS2(s:float,alpha:float)->float:
    k,c,b = .504,.14601,4.
    p     = s**2

    F1x = 1 + k - k/(1+mu*p/k)
    F0x = 1 + k - k/(1+(mu*p+c)/k)

    f = (1-alpha**2) / (1 + alpha**3 + b*alpha**6)

    return F1x + f*(F0x-F1x)


##############################################################################
# Helper Functions
#----------------

def leg(n:int,x:float)->float:
    """
    Value of nth legendre polynomial evaluated at x
    """
    v = np.zeros(n+1)
    v[n] = 1 # n = 6 ===> [0,0,0,0,0,1]
    return np.polynomial.legendre.Legendre(v)(x)

def transformed_s(s:float)->float:
    """
    Properties: s=0 -> -1 , s=inf -> 2

    Eq # 2 in mBEEF paper
    """
    kappa = 0.804
    mu    = 10./81
    q     = kappa / mu # Pade approximant to PBEsol Fx(s)
    return 2*s**2 / (q+s**2) - 1

def transformed_alpha(alpha:float)->float:
    """
    Properties: a=0 -> 1, a = 1 -> 0, a = inf -> -1

    Eq # 3 in mBEEF paper

    MULTIPLIED BY NEGATIVE ONE!!!
    """
    return -(1-alpha**2)**3/(1+alpha**3 + alpha**6)


def LegendreProduct(s:float,alpha:float,M:int,N:int)->float:
    """
    s     - reduced density gradient
    alpha - reduced kinetic energy density
    M,N   - L. polynomials

    Eq # 4 in mBEEF paper
    """
    return leg(M,transformed_s(s)) * leg(N,transformed_alpha(alpha))


def Enhancement(s:float,alpha:float,A:array)->float:
    """
    A - an MxN matrix with rows corresponding to s basis functions coefficients
        and columns corresponding to alpha basis function coefficients

    Eq # 5 in mBEEF paper
    """
    N,M = A.shape
    P = array([[LegendreProduct(s,alpha,j,i) for j in range(N)] for i in range(M)])
    return np.sum(np.multiply(A,P))

########
# Parsing
#-------
def readcsv(fi:str,x_size:int=8) -> csv_out_type:
    """
    Extract element,n_unit,target,Eatom,Ebulk  from CSV

    Inputs
        fi     - path to CSV file
        x_size - dimensions of exchange contribution matrix being parsed
    Outputs
        N x (m^2) matrices:
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
    XAtom,XBulk = [np.empty((0,x_size**2))]*2

    # Parse
    with open(fi, 'r') as f:
        reader = DictReader(f, delimiter=',', quotechar='"')
        for row in reader:

            # add the new row's values to the correct list
            for c in simpleCols:
                simple[c].append(row[c])

            # 'flatten' the mxm matrices into vectors, then 'append'
            xAtom,xBulk = [array(literal_eval(row[y])).reshape((1,x_size**2))
                                for y in ['xAtom','xBulk']]

            XAtom = np.vstack((XAtom,xAtom))
            XBulk = np.vstack((XBulk,xBulk))

    # Postrocess the resulting lists
    element,n_unit,target,Eatom,Ebulk = simple.values()

    return (element
           ,array(n_unit).astype(int)
           ,array(target).astype(float)
           ,array(Eatom).astype(float)
           ,array(Ebulk).astype(float)
           ,XAtom,XBulk)

###############
# Visualization
#--------------
def plot(f      : binary = add
        ,label  : str                           = ''
        ,color  : str                           = 'k'
        ,mgga   : bool                          = False
        ) -> None:

    ss     = np.arange(0.,3.5,0.05)
    alphas = np.array([0,1]) if mgga else [1]

    styles = ['-',':','--','-.']

    for style,alpha in zip(styles,alphas):
        ys = [f(s,alpha) for s in ss]
        lab = label + r' ($\alpha$=%d)'%alpha if mgga else label
        plt.plot(ss,ys,color=color,linestyle=style,label=lab,alpha=1 if mgga else 0.3)

# Constraints

class LinConstraint(object):
    """
    Force the functional to take a value at a particular s,alpha point
    """

    weight = 100

    def __init__(self,s:float,alpha:float,val:float)->None:
        self.s      = s
        self.alpha  = alpha
        self.val    = val

    @staticmethod
    def add(X:array,t:array,lst:List['LinConstraint'])->Tuple[array,array]:
        """
        Update an X matrix and target vector with a list of Constraints
        """
        m  = int(X.shape[1]**(0.5))
        print('m = %d'%m)
        for const in lst:
            cx = array([LegendreProduct(const.s,const.alpha,j,i)
                                    for i in range(m) for j in range(m)])
            X = np.vstack((X,cx * const.weight))
            t = np.concatenate((t,[const.val * const.weight]))
        return X,t

# Constraint Definitions
#-----------------------

cLDA    = [LinConstraint(0,1,1)]
cLiebOx = [LinConstraint(1e10,1,1.804),LinConstraint(1e10,1e10,1.804)]
cGD     = [LinConstraint(1e10,0,1.174)]#,LinConstraint(1e10,1e10,1.174)]
cMS2    = LinConstraint(1e10,2,1.5)


# # ,[LegendreProduct(1e10,2.,j,i) for i in range(m) for j in range(m)] -> 1.5
# # ,[LegendreProduct(1e10,0,j,i) for i in range(m) for j in range(m)]]) -> 1.4



###############################################################################
def main(pth    : str
        ,m      : int   = 5
        ,x_size : int   = 8
        ,reg    : float = 0.001
        ) -> None:
    """
    Inputs:
    pth    - filepath to CSV file containing cols with elements and xc contribs
    m      - number of Legendre basis functions (for both s and alpha)
    x_size - length of the square matrix of Legendre coefficients in CSV file
    reg    - regularization for least square fitting
    weight - weight of contraints in least square fitting
    """
    assert m <= x_size # max is 8x8

    # Parse file
    #------------
    elems,natoms,Ece,Eatom,Ebulk,Xatom,Xbulk = readcsv(pth,x_size)

    # Get real target
    #----------------
    Xbulk_norm = Xbulk / natoms[:,None]
    Ebulk_norm = Ebulk / natoms
    Exatom = np.dot(Xatom,beefcoeff)
    Exbulk = np.dot(Xbulk_norm,beefcoeff)
    dE     = (Ebulk_norm - Exbulk) - (Eatom - Exatom) # Non exchange contribution of cohesive energy (calculated by DFT, assumed valid)
    target = Ece - dE          # 'True' difference E_x (bulk - atom)

    # Get real fitting inputs
    #------------------------
    remove = filter(lambda x: x % x_size >= m or x >= x_size*m,range(x_size**2))
    X      = np.delete(Xbulk_norm - Xatom,list(remove),1)

    # Apply Constraints
    #------------------
    constraints = cLDA+cGD+cLiebOx
    X,target = LinConstraint.add(X,target,constraints)

    #Least-squares fit (rcond controls over/under fitting)
    #----------------------------------------------------
    fit,resid,rank,singvals = np.linalg.lstsq(X, target, rcond=reg)

    # Calculate residual
    #-------------------
    res = np.dot(X, fit) - target
    print('target %s '%target)
    print('residual %s'%res)
    #Effective number of degrees of freedom gives a scale of the ensemble
    #temperature T something to play with, controls the Spaghetti spread
    #----------------------------------------------------------------------
    Ndeg = 2.5                          # Equation 8 (for actual formula)
    min_cost = 0.5 * np.sum(res**2)     # Equation 9
    T   = 2. * min_cost / Ndeg          # Equation 11b
    A   = np.dot(X.T,X) / T
    l,U = np.linalg.eigh(A)
    L   = np.diag(1./np.sqrt(l))
    M   = np.dot(U,L)  # ensemble information

    # Make function from fit
    #-----------------------
    f = lambda s,a: Enhancement(s,a#beefcoeff.reshape(8,8))
                        ,fit.reshape((m,m)))

    # Make plots
    #-----------
    plot(fxPBE,'PBE','g')
    plot(fxRPBE,'RPBE','r')
    #plot(fxPBEsol,'PBEsol','b')
    plot(fxSCAN,'SCAN','k',mgga=True)
    plot(fxMS2,'MS2','purple',   mgga=True)
    plot(f,'Cohesive Energy Fit','b',mgga=True)

    # Manager overall plot settings
    #------------------------------
    plt.xlabel('Reduced density gradient')
    plt.ylabel('Exchange enhancement')
    plt.legend(loc='best')
    #plt.ylim((0.6,2.0))


if __name__=='__main__':
    # Validate inputs
    assert exists(argv[1])

    subplot = False

    LinConstraint.weight = 10000

    if subplot:
        # Subplot stuff
        vals = [0,0.01,1]
        for i in range(len(vals)):
            plt.subplot(len(vals),1,i+1)
            main(argv[1],3,reg=vals[i])
            plt.title('regularization = %s'%vals[i])
        plt.subplots_adjust(hspace=0.6)
    else:
        # just one plot
        main(argv[1],2,reg=0.1)
    plt.show()
