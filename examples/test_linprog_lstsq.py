import numpy as np # type: ignore
from functionals.fit.utilities import lst_sq_linprog
"""
Minimize Ax = b using linear programming

A = [2 1
     0 1]

b = [3 1]

Solution should be [1,1]
"""



if __name__=='__main__':
    A = np.array([[2,1],[0,1]])
    b = np.array([3,1])

    # But now, we can add bounds to our x variables
    bounds =  [(-10.,10.),(-10,10)]

    # And we can have linear eq/ineq constraints on the solution
    A_eq = np.array([[2,2],[4,4]]) # (things that happen to be true
    b_eq = np.array([4,8.001])     # of the unique solution)

    A_ub = np.array([[2,2],[4,4]]) # (trivially true, within tolerance,
    b_ub = np.array([3.99,9])      # given the EQ constraints)

    print('\nexpecting [1. 1.]')
    print(lst_sq_linprog(A,b,bounds,A_eq,b_eq,A_ub,b_ub))
