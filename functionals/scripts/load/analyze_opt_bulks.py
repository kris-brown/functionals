from typing import Tuple as T, Optional as O
from json import loads, dumps
import numpy as np


def analyze_opt_bulks(pth: str, vols_: str, engs_: str, struct: str
                      ) -> T[O[float], O[float], O[str], O[float], O[float]]:
    '''Analysis on serialized info about 5 jobs centered on optimum.'''
    import ase
    from ase.eos import EquationOfState
    # DEBUG = False
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    zipped = zip(loads(vols_), loads(engs_))
    s_vols, s_engs = zip(*sorted(zipped))

    assert len(s_vols) == 5
    dxs = set(np.round(np.diff(s_vols), 7))
    if len(dxs) > 1:
        print('Bad spacings ', pth, dxs, vols_)
        return None, None, None, None, None
    dx = dxs.pop()
    stencil = (-1 * s_engs[4] + 16 * s_engs[3] - 30 * s_engs[2] +
               16 * s_engs[1] - s_engs[0]) / (12 * dx**2)

    bulkmod = stencil * s_vols[2] * ev_a3_to_gpa  # APPROXIMATE
    # Quadratic fit
    fit = np.polyfit(s_vols, s_engs, 2)
    lstsq_abc = dumps(fit.tolist())
    bm_lstsq = float(-fit[1]) * ev_a3_to_gpa  # APPROXIMATE

    eos = EquationOfState(s_vols, s_engs)
    try:
        _, __, eosbm = eos.fit()  # type: T[float,float,float]
    except Exception as e:
        print("%s optbulkeos has error %s" % (pth, e))
        return None, bm_lstsq, lstsq_abc, None, None
    eosbm = eosbm / ase.units.kJ * 1.0e24

    vol = s_vols[2] * (4 if 'prim' in struct else 1)
    return bulkmod, bm_lstsq, lstsq_abc, eosbm, vol
