import unittest
import functionals.fit.math as math
import numpy as np


class TestFunc(unittest.TestCase):
    '''
    Is there a way to test whether h_norm_const is correct?
    '''

    def test_chain_rule(self) -> None:
        '''Simple example with known solution.'''
        def f(x): return x**2 + np.exp(x)
        def f_(x): return 2*x + np.exp(x)
        def f__(x): return 2 + np.exp(x)
        def f___(x): return np.exp(x)
        def g(x): return 3*x**3 + 1
        def g_(x): return 9*x**2
        def g__(x): return 18*x
        def g___(x): return 18

        # d2(f(g(x)))
        dfgx2_dx2 = math.chain_rule2(f__, f_, g__, g_, g)(0)
        self.assertEqual(dfgx2_dx2, 0)
        dfgx2_dx2 = math.chain_rule2(f__, f_, g__, g_, g)(0.1)
        self.assertAlmostEqual(dfgx2_dx2, 8.55669, places=4)
        dfgx3_dx3 = math.chain_rule3(f___, f__, f_, g___, g__, g_, g)(0.1)
        self.assertAlmostEqual(dfgx3_dx3, 87.4831, places=4)


if __name__ == "__main__":
    unittest.main()
