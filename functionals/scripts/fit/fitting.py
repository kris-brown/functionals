# External
from typing import Optional as O, Tuple as T
# Internal
from functionals.fit           import Fit

###############################################################################
def fitting(data            : str,
            constraints     : str,
            nlconstraints   : str,
            maxit           : int,
            n               : int,
            bound           : float,
            ifit            : bool,
            gridden         : int,
            bm_weight       : float,
            lat_weight      : float
           ) -> T[O[str],float,O[int],O[str],str]:

    return Fit(data=data, constraints=constraints,
                nlconstraints=nlconstraints,maxit=maxit,n=n,bound=float(bound),
                ifit=ifit,gridden=gridden,bm_weight=float(bm_weight),
                lat_weight=float(lat_weight)).fit_result()
