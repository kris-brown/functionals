# External
from typing      import Tuple,List,Dict
from sys         import argv
from csv         import DictReader,reader
from collections import defaultdict
from ast         import literal_eval
from os.path     import exists
from matplotlib  import rc             # type: ignore
from numpy       import array          # type: ignore
from matplotlib  import pyplot as plt  # type: ignore
import numpy as np                     # type: ignore
# Internal Modules
from functionals.functional import PBE,RPBE,PBEsol,SCAN,MS2,BEEF,FromMatrix,beefcoeff
from functionals.utilities  import LegendreProduct
###############################################################################
plt.rcParams.update({"pgf.texsystem": "pdflatex"})
np.set_printoptions(precision=2,linewidth=120,floatmode='fixed',suppress=True)
csv_out_type = Tuple[List[str],array,array,array,array,array,array]

##############################################################################


########
# Parsing
#-------
def readcsv(fi:str,x_size:int=8) -> csv_out_type:
    """
    Extract element,n_unit,target,Eatom,Ebulk from CSV

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


#[LegendreProduct(1e10,2.,j,i) for i in range(m) for j in range(m)] -> 1.5
#[LegendreProduct(1e10,0,j,i) for i in range(m) for j in range(m)]]) -> 1.4



###############################################################################
def main(pth    : str,
         m      : int   = 5,
         x_size : int   = 8,
         reg    : float = 0.001
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
    Exatom     = np.dot(Xatom,beefcoeff)
    Exbulk     = np.dot(Xbulk_norm,beefcoeff)
    dE         = (Ebulk_norm - Exbulk) - (Eatom - Exatom) # Non exchange contribution of cohesive energy (calculated by DFT, assumed valid)
    target     = Ece - dE                                 # 'True' difference E_x (bulk - atom)

    # Get real fitting inputs
    #------------------------
    remove = filter(lambda x: x % x_size >= m or x >= x_size*m,range(x_size**2))
    X      = np.delete(Xbulk_norm - Xatom,list(remove),1)

    # Apply Constraints
    #------------------
    constraints = cLDA + cGD + cLiebOx
    X,target    = LinConstraint.add(X,target,constraints)

    #Least-squares fit (rcond controls over/under fitting)
    #----------------------------------------------------
    fit,resid,rank,singvals = np.linalg.lstsq(X, target, rcond=reg)

    # Calculate residual
    #-------------------
    res = np.dot(X, fit) - target
    print('target %s ' %target)
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
    custom = FromMatrix(fit.reshape((m,m)),name='Cohesive Energy Fit')
    # Make plots
    #-----------
    ax = plt.gca()
    PBE.plot(ax,'g')
    RPBE.plot(ax,'r')
    SCAN.plot(ax,'k')
    MS2.plot(ax,'purple')
    BEEF.plot(ax,'orange')
    custom.plot(ax,'b')

    # Manager overall plot settings
    #------------------------------
    plt.xlabel('Reduced density gradient')
    plt.ylabel('Exchange enhancement')
    plt.legend(loc='best')
    #plt.ylim((0.6,2.0))
    plt.show()
    plt.savefig('spaghetti.pdf',bbox_inches='tight')


if __name__=='__main__':
    # Validate inputs
    #----------------
    assert exists(argv[1])

    # Parameters
    #-----------
    LinConstraint.weight = 10000

    # Main
    #-----
    main(argv[1],2,reg=0.1)


# if subplot:
#     # Subplot stuff
#     vals = [0,0.01,1]
#     for i in range(len(vals)):
#         plt.subplot(len(vals),1,i+1)
#         main(argv[1],3,reg=vals[i])
#         plt.title('regularization = %s'%vals[i])
#     plt.subplots_adjust(hspace=0.6)
