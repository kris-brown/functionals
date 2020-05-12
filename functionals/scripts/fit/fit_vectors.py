import pdb
from typing import Tuple as T, Optional as O
import numpy as np


def fit_vectors(engs_: str, vols_: str, contribs_: str, coefs_: str,
                volume: float, name: str
                ) -> T[O[str], O[str]]:
    '''Ax=b fit vectors for experimental CE, BM, and volume.

    Quadratic approximation around center point: Eng = Av^2 + Bv + C

    E[i](x) = Enonxc[i] + C[i]@x
    Enonxc[i] = E[i] - C[i]-CalcCoefs

    5 point appx for 2nd derivative (i.e. E''= 2A):
        -E[+2] + 16E[+1] - 30E[0] + 16E[-1] - E[-2] / 12h^2

    4 point appx for 1st derivative (i.e. E'=2Av + B.... or B=E'-2Av):
        -E[+2] + 8E[+1] - 8E[-1] + E[-2] / 12h

    Bulk modulus is equal to A = E''/2

    optimal volume: v_opt = -B/2A
                    v_opt = -(E'-2Av)/E''
                    v_opt = v - E'/E''

    - [[eV/A³ to GPa]] = (10**-9 GJ/J) * (1.602 * 10**-19 J/eV) * (10**10 A/m)**3
    - GPa = GJ/m³
    '''
    from json import dumps, loads
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    engs, vols, contribs, coefs = [np.array(loads(x)) for x in [
        engs_, vols_, contribs_, coefs_]]
    if len(engs) != 5:
        return None, None
    dx = np.mean(np.diff(vols))
    assert np.max(np.diff(vols)-dx) < 0.0001, vols  # check even spacing
    volume = float(volume)

    class AB(object):
        '''A pair A,b such that Ax=b for some vector of coefficients.'''

        def __init__(self, A: np.array = None, b: float = None) -> None:
            self.A = A if A is not None else np.zeros(64,)
            self.b = b or 0.

        def __str__(self) -> str: return dumps([self.A.tolist(), self.b])
        def __add__(self, ab: 'AB') -> 'AB': return AB(self.A +
                                                       ab.A, self.b+ab.b)
        def __sub__(self, ab: 'AB') -> 'AB': return AB(self.A -
                                                       ab.A, self.b-ab.b)

        def __rmul__(self, z: float) -> 'AB': return AB(self.A*z, self.b*z)
        def __truediv__(self, z: float) -> 'AB': return AB(self.A/z, self.b/z)
        def __matmul__(
            self, v: np.array) -> float: return float(self.A @ v + self.b)

    Enonxc = [e-contrib@coefs for e, contrib in zip(engs, contribs)]
    E = [AB(contrib, enonxc) for contrib, enonxc in zip(contribs, Enonxc)]
    ddE = (-1*E[4] + 16*E[3] - 30*E[2] + 16*E[1] - E[0])/(12*dx**2)
    Bulkmod = float(volume) * ev_a3_to_gpa * ddE

    # dE = (-1*E[4] + 8*E[3] - 8*E[1] + E[0])/(12*dx)
    # old_Vopt = AB(b=float(volume)) - dE/(ddE@coefs)

    JohannesMagicQ = -1/(2*(ddE@coefs)*(dx**2))
    v1prefactor = JohannesMagicQ*(-2*volume - dx)
    v2prefactor = JohannesMagicQ*4*volume
    v3prefactor = JohannesMagicQ*(-2*volume + dx)

    Vopt = v1prefactor * E[1] + v2prefactor * E[2] + v3prefactor * E[3]
    if name == 'Nb':
        import pdb
        pdb.set_trace()
    return str(Bulkmod), str(Vopt)
