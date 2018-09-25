# External
from typing      import List
from sys         import argv
from os.path     import exists
from matplotlib  import pyplot as plt  # type: ignore
import numpy as np                     # type: ignore

# Internal Modules
from functionals.constraint import Constraint,cLiebOx,cLDA
from functionals.data       import CohesiveData,Data
from functionals.fit        import Fit
from functionals.functional import PBE, RPBE, SCAN, MS2, BEEF
###############################################################################
# Config
#--------
plt.rcParams.update({"pgf.texsystem": "pdflatex"})
np.set_printoptions(precision=2, linewidth=120, floatmode='fixed', suppress=True)

##############################################################################
def main(pth: str, m: int = 8, verbose: bool  = False) -> None:
    """
    Inputs:
    pth    - filepath to CSV file containing cols with elements and xc contribs
    m      - number of Legendre basis functions (for both s and alpha)
    """
    data = [CohesiveData(pth)] # type: List[Data]
    constraints = [cLiebOx,cLDA] # type: List[Constraint] ### [cLiebOx,cLDA]
    fit = Fit(m,data,constraints).functional()


    # Make plots
    #-----------
    ax = plt.gca()
    PBE.plot(ax,'g')
    RPBE.plot(ax,'r')
    SCAN.plot(ax,'k')
    MS2.plot(ax,'purple')
    BEEF.plot(ax,'orange')
    fit.plot(ax,'blue')
    # Manager overall plot settings
    #------------------------------
    plt.xlabel('Reduced density gradient')
    plt.ylabel('Exchange enhancement')
    plt.legend(loc='best')
    plt.show()
    plt.savefig('spaghetti.pdf',bbox_inches='tight')

if __name__=='__main__':

    # Validate inputs
    #----------------
    assert exists(argv[1]), "%s doesn't point to a CSV file "%argv[1]

    # Parameters
    #-----------
    verbose  = True

    # Main
    #-----
    main(argv[1],m = 3,verbose=verbose)
