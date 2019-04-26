from json import load, dumps, dump
from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint # type: ignore
from numpy import array,vstack,concatenate,zeros # type: ignore
from numpy.linalg   import lstsq                            # type: ignore

from io         import StringIO
from contextlib import redirect_stdout
from csv        import reader

'''
Just a file to play around with, nothing serious

'''
# Constants
###########
params  = dict(  verbose                  = 3,
                 disp                     = True,
                 maxiter                  = 100000,
                 xtol                     = 1e-10,#10,
                 barrier_tol              = 1e-4,#10,
                 gtol                     = 1e-3,#10,
                 initial_tr_radius        = 1.0,
                 initial_constr_penalty   = 1.0,
                 initial_barrier_parameter= 0.1)

def main() -> None:

    # Extract linear constraints
    #############################

    c_A  = [[3,4],[2,-4]]
    c_lb = [1,-1]
    c_ub = [2,2]

    cons   = LinearConstraint(c_A, c_lb, c_ub, keep_feasible = False)

    # Extract nonlinear constraints
    ###############################
    all_cons = [cons]

    # Extract data to be fit
    ########################

    initguess = [10,10]

    # Define cost functions
    #######################

    def f(x : array) -> float:
        """Objective function for optimization: sum of squares of residual"""
        return x @ x

    # Call SciPy Minimize
    #####################

    bounds = Bounds([0.,-2.], [2.,2.] )

    redir = StringIO()

    with redirect_stdout(redir):

        res = minimize(f, initguess, method = 'trust-constr',
                        constraints = all_cons,jac='2-point',hess=lambda x: [0,0],
                        bounds = bounds, options = params)

    # Serialize results
    ###################
    dt = int(res.execution_time);
    x  = dumps(res.x.tolist())

    lines  = [line[1:-1] for line in redir.getvalue().split('\n')[2:-4]
                if line.count('|')==11]
    output = [(float(n),float(o),float(c))
                for n,_,_,o,_,_,c,_,_,_ in reader(lines,delimiter='|')]

if __name__=='__main__':
    main()
