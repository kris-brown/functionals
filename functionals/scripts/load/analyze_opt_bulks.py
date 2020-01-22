from typing import Tuple as T
from json import loads, dumps
import numpy as np


def analyze_opt_bulks(pth: str, vols_: str, engs_: str) -> T[float, float, str, float]:
    '''Analysis on serialized info about 5 jobs centered on optimum.'''
    import ase
    from ase.eos import EquationOfState
    # DEBUG = False
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    serial = [vols_, engs_]
    s_vols, s_engs = map(loads, serial)

    assert len(s_vols) == 5
    dxs = set(np.round(np.diff(s_vols), 7))
    if len(dxs) > 1:
        print('Bad spacings ', pth, dxs, vols_)
        assert False, (pth, dxs, vols_)
    dx = dxs.pop()
    stencil = (-1*s_engs[4] + 16*s_engs[3] - 30*s_engs[2]
               + 16*s_engs[1] - s_engs[0])/(12*dx**2)

    bulkmod = stencil * s_vols[2] * ev_a3_to_gpa  # APPROXIMATE
    # Quadratic fit
    fit = np.polyfit(s_vols, s_engs, 2)
    lstsq_abc = dumps(fit.tolist())
    bm_lstsq = float(-fit[1]) * ev_a3_to_gpa  # APPROXIMATE

    eos = EquationOfState(s_vols, s_engs)
    _, __, eosbm = eos.fit()  # type: T[float,float,float]
    eosbm = eosbm / ase.units.kJ * 1.0e24

    if s_vols[0] == 43.420321409573106:
        import pdb
        pdb.set_trace()
    return bulkmod, bm_lstsq, lstsq_abc, eosbm
