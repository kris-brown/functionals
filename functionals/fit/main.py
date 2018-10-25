# External
from typing      import List,Tuple
from sys         import argv
from os.path     import exists

from matplotlib  import pyplot as plt  # type: ignore
import numpy as np                     # type: ignore

# Internal Modules
from functionals.fit.constraint import Constraint,FxConstraint,cLiebOx,cLDA,cSCAN11,cPos
from functionals.fit.data       import CohesiveData,Data
from functionals.fit.fit        import LinFit,NonLinFit
from functionals.fit.functional import RPBE, SCAN, MS2, BEEF
###############################################################################
# Config
#--------
plt.rcParams.update({"pgf.texsystem": "pdflatex"})
np.set_printoptions(precision=2, linewidth=120, floatmode='fixed', suppress=True)
ax = plt.gca()
################################################################################
def main(pth        : str,
         basis      : int                = 8,
         nonlin     : bool               = True,
         norm       : float              = 0.1,
         initfit    : bool               = True,
         bound      : float              = 0.1,
         verbose    : bool               = False,
         maxiter    : int                = 1000,
         constgrid  : Tuple[int,int,int] = (0,4,10)
         ) -> None:
    """
    Inputs:
    pth    - filepath to CSV file containing cohesive energy data
    basis  - number of Legendre basis functions (for both s and alpha)
    nonlin - use nonlinear fitting
    norm   - regularization term (only used if nonlin)
    initfit- how to initialize nonlinear fitting (only used if nonlin)
    bound  - for all non-constant terms of fitted coefs, bound the magnitude
    """
    data        = CohesiveData(pth)
    constraints = [cLiebOx,cLDA,cSCAN11,cPos] # type: List[Constraint]

    for c in constraints:
        if isinstance(c,FxConstraint):
            c.ss = np.linspace(*constgrid);

    fitter = NonLinFit if nonlin else LinFit
    args   = {'norm':norm,'initfit':initfit} if nonlin else {}
    fit    = fitter(basis,[data],constraints,bound=bound,**args)
    fx     = fit.functional(maxiter=maxiter,verbose=verbose)

    # Make plots
    #-----------
    colors = ['r', 'k', 'purple', 'g', 'blue']
    fxs    = [RPBE, SCAN, MS2, BEEF, fx]

    for fx,col in zip(fxs, colors):
        fx.plot(ax, col)

    loss = fx.resid(data)[0] - BEEF.resid(data)[0]


    # Manager overall plot settings
    #------------------------------
    xmin, xmax, ymin, ymax = ax.axis()


    fmtargs = [basis,'non' if nonlin else '',
               ',norm %s, ifit %s'%(norm,initfit) if nonlin else '',
               bound,fit.runtime]

    plt.xlabel('Reduced density gradient')
    plt.ylabel('Exchange enhancement')
    plt.text(xmin,ymax,'"Loss" rel. to BEEF of %d'%int(loss),verticalalignment='top')
    plt.title('Basis {},{}lin{},bound={},time={} '.format(*fmtargs))
    plt.legend(loc = 'best')
    plt.show()
    plt.savefig('spaghetti.pdf',bbox_inches = 'tight')

if __name__=='__main__':

    # Validate inputs
    #----------------
    assert len(argv)==2, 'Need to provide path to CSV file as only CLI argument'
    assert exists(argv[1]), "%s doesn't point to a CSV file "%argv[1]

    # Parameters
    #-----------
    verbose   = False
    nonlin    = True
    norm      = 0.001
    initfit   = False
    basis     = 8
    constgrid = (0,5,10)
    bound     = 0.1 # max |value| of any element in the fitted matrix, except for 1st element (the offset)
    maxiter   = 2000
    # Main
    #-----
    main(pth = argv[1],basis = basis, nonlin = nonlin, constgrid=constgrid,
         norm = norm, bound = bound, initfit = initfit, maxiter = maxiter,
         verbose = verbose)
