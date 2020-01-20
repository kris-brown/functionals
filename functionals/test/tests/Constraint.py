import unittest
import functionals.fit.constraint as con
import numpy as np

class TestFunc(unittest.TestCase):
    '''
    Is there a way to test whether h_norm_const is correct?
    '''

    def test_chain_rule(self) -> None:
        '''Simple example with known solution.'''
        f = lambda x: x**2 + np.exp(x)
        f_ = lambda x: 2*x + np.exp(x)
        f__ = lambda x: 2 + np.exp(x)
        f___ = lambda x: np.exp(x)
        g = lambda x: 3*x**3 + 1
        g_ = lambda x: 9*x**2
        g__ = lambda x: 18*x
        g___ = lambda x: 18
        dfgx2_dx2 = con.chain_rule2(f__,f_,g__,g_,g,0)
        self.assertEqual(dfgx2_dx2,0)
        dfgx2_dx2 = con.chain_rule2(f__,f_,g__,g_,g,0.1)
        self.assertAlmostEqual(dfgx2_dx2,8.55669,places=4)
        dfgx3_dx3 = con.chain_rule3(f___,f__,f_,g___,g__,g_,g,0.1)
        self.assertAlmostEqual(dfgx3_dx3,87.4831,places=4)


if __name__ == "__main__":
    unittest.main()
