from typing                    import Tuple as T
from json                      import load, dumps, dump, loads
from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint # type: ignore
from numpy                     import array,vstack,multiply,concatenate as concat,zeros,logspace, empty, exp,inf # type: ignore
from numpy.linalg              import lstsq,norm                            # type: ignore
from numpy.polynomial.legendre import Legendre                              # type: ignore
from io         import StringIO
from contextlib import redirect_stdout
from csv        import reader
################################################################################

# Constants
###########
params  = dict(  verbose                  = 3,
                 disp                     = True,
                 maxiter                  = 100000,
                 xtol                     = 1e-10,
                 barrier_tol              = 1e-4,
                 gtol                     = 1e-3,
                 initial_tr_radius        = 1.0,
                 initial_constr_penalty   = 1.0,
                 initial_barrier_parameter= 0.1)
# Helpers
###########
def flatten(lol : list) -> list: return [item for sublist in lol for item in sublist]

def leg(n : int,x : float) -> float:
    """ Value of nth legendre polynomial evaluated at x """
    v = zeros(n+1)
    v[n] = 1 # n = 6 ===> [0,0,0,0,0,1]
    return Legendre(v)(x)

def transformed_s(s : float) -> float:
    """ Properties: s=0  ->  -1 , s=inf  ->  1 (Eq # 2 in mBEEF paper) """
    assert s >= 0
    kappa = 0.804
    mu    = 10./81
    q     = kappa / mu # Pade approximant to PBEsol Fx(s)
    return 2 * s**2 / (q+s**2) - 1

def transformed_alpha(alpha:float, msb : float) -> float:
    """
    Properties: a=0 -> 1, a = 1 -> 0, a = inf -> -1
    Eq # 3 in mBEEF paper
    MULTIPLIED BY NEGATIVE ONE!!!
    """
    assert alpha >= 0
    return -(1-alpha**2)**3/(1+alpha**3 + msb*alpha**6)

def LegendreProduct(s      : float,
                    alpha  : float,
                    a1     : float,
                    msb    : float,
                    M      : int,
                    N      : int
                   ) -> float:
    """
    s     - reduced density gradient
    alpha - reduced kinetic energy density
    M,N   - L. polynomials
    Eq # 4 in mBEEF paper
    """
    s   = max(s,1e-30); alpha = max(alpha,1e-30)
    gx  = 1 - exp(-a1*s**(-0.5))
    return gx * leg(M,transformed_s(s)) * leg(N,transformed_alpha(alpha,msb))

def fx(s:float,alpha:float,A:array,a1:float,msb:float) -> float:
    return sum([A.reshape((8,8))[i,j] * LegendreProduct(s,alpha,a1,msb,i,j)
                    for i in range(8) for j in range(8)])

def linconst(const : dict, a1 : float, msb : float, const_den : int) -> T[bool,array,array]:
    ''' Turns a linear constraint into an Ax = b matrix '''
    w = const['weight']
    assert w > 0
    grid = logspace(-1,0.8,const_den).tolist()

    val = float(const['val'])

    if   const['kind'] == 'eq':  eq = True  ; sign = 1
    elif const['kind'] == 'lt':  eq = False ; sign = 1
    elif const['kind'] == 'gt':  eq = False ; sign = -1
    else: raise ValueError

    if const['vec']:
        c_A = loads(const['vec'])
        vals = [val]
    else:
        s_, a_ = const['s'],const['alpha']
        srange = [float(s_)] if s_ is not None else grid
        arange = [float(a_)] if a_ is not None else grid
        card   = len(srange) * len(arange)
        #w      = w / card # reduce weight if many data points
        vals   = [val for _ in range(card)]
        c_A    = []
        for s in srange:
            for a in arange:
                c_A.append(flatten([[LegendreProduct(s,a,a1,msb,i,j)
                                        for j in range(8)]
                                            for i in range(8)]))
    return eq, w * array(c_A)  * sign, w * array(vals)  * sign

def main() -> None:

    with open('metadata.json','r') as fi: params = load(fi)['params']

    # Extract linear constraints
    #############################
    c_A_eq,c_A_lt,c_b_eq,c_b_lt = empty((0,64)),empty((0,64)),empty((0,)),empty((0,)),

    with open('constraints.json','r') as fi:
        for const in load(fi):
            eq,A,b = linconst(const,params['a1'],params['msb'],params['cd'])
            if eq:
                c_A_eq = vstack((c_A_eq,A))
                c_b_eq = concat((c_b_eq,b))
            else:
                c_A_lt = vstack((c_A_lt,A))
                c_b_lt = concat((c_b_lt,b))

            neginf = array([-inf for _ in c_b_lt])

            c_A = vstack((c_A_eq,c_A_lt))
            c_lb = concat((c_b_eq,neginf))
            c_ub = concat((c_b_eq,c_b_lt))
            cons = LinearConstraint(c_A, c_lb, c_ub, keep_feasible = False)

    # Extract nonlinear constraints
    ###############################


    # Extract data to be fit
    ########################

    with open('data.json','r') as fi:
        ceX,bmX,lX,ceY,bmY,lY = map(array,load(fi))

    X = vstack((ceX,bmX,lX))
    Y = concat((ceY,bmY,lY))

    initguess = lstsq(X,Y,rcond=None)[0]

    # Define cost functions
    #######################

    def f(x : array) -> float:
        """Objective function for optimization: sum of squares of residual"""
        resid = X @ x - Y
        return resid @ resid

    def df(x : array) -> array:
        """Jacobian of objective function (derived below in Einstein notation)"""
        residual = X @ x - Y
        return 2 * X.T @ residual

    def hess(x:array,*args:array) -> array:
        return zeros((64,64))

    # Call SciPy Minimize
    #####################

    bounds = Bounds([0.] + [-2.] * 63, [2.] + [2.] * 63)

    redir = StringIO()

    with redirect_stdout(redir):
        res = minimize(f, initguess, method = 'trust-constr', jac = df,
                        constraints = [cons], hess = hess,
                        bounds = bounds, options = params)

    # Serialize results
    ###################
    dt = int(res.execution_time);
    x  = dumps(res.x.tolist())

    # Parse output csv
    ###################
    # Header: niter |f evals|CG iter|  obj func   |tr radius |   opt    |  c viol  | penalty  |barrier param|CG stop|
    lines  = [line[1:-1] for line in redir.getvalue().split('\n')[2:-4]
                if line.count('|') == 11]

    output = [(int(n),float(o),float(c))
                for n,_,_,o,_,_,c,_,_,_ in reader(lines,delimiter='|')]

    ;breakpoint()

    with open('output.json','w') as fi:
        dump(output,fi)

    with open('result.json','w') as fi:
        dump([dt,x],fi)

if __name__=='__main__':
    main()


"""
BONUS
-----

Derivation of Jacobian in Einstein notation:

        Capital letters = vectors, matrices in [brackets]

        Fi  = residual                                                      --- (DIM: m x 1)
        Aij = basis function coefficients (cols) for each data point (row)  --- (DIM: m x n)
        Rj  = vector of weights corresponding to basis functions            --- (DIM: n x 1)
        Yi  = vector of targets for dot product of rows in Aij and Rj       --- (DIM: m x 1)
        δij = Kroenecker delta

        let: Fj = [A]ji * Ri - Yj
        d(Fj)/d(Rk) = [A]ji * d(Ri)/d(Rk) = [A]ji * δik = [A]jk

        d(FjFj)/d(Rk) = 2 * d(Fj)/d(Rk) * Fj
                      = 2 * [A]jk * Fj

Derivation of constant hessian:
        To see this, take derivative of jac result w/r/t some Rm:

        d(FjFj)/d(Rk)d(Rm) = 2 * [A]jk * d(Fj)/d(Rm) =  2 * [A]jk * [A]jm
"""
