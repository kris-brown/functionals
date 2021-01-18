from typing import Tuple as T, Optional as O
from json import loads, dumps
import numpy as np


def analyze_opt_bulks(pth: str, vols_: str, engs_: str, struct: str
                      ) -> T[O[float], O[float], O[float], O[str]]:
    '''Analysis on serialized info about 5 jobs centered on optimum.'''

    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    zipped = zip(loads(vols_), loads(engs_))
    s_vols, s_engs = zip(*sorted(zipped))

    assert len(s_vols) == 5

    if len(set(np.round(np.diff(s_vols), 3))) > 1:
        print('Bad spacings ', pth, vols_)
        return None, None, None, None

    dx = np.mean(np.diff(s_vols))
    stencil = (-1 * s_engs[4] + 16 * s_engs[3] - 30 * s_engs[2] +
               16 * s_engs[1] - s_engs[0]) / (12 * dx**2)

    bulkmod = stencil * s_vols[2] * ev_a3_to_gpa  # APPROXIMATE
    if bulkmod < 0:
        print(pth, ' -- bulkmod below 0')
        return None, None, None, None

    Q = -1 / (2 * stencil * dx**2)
    v1pre = Q * (s_vols[2] * -2 - dx)
    v2pre = Q * 4 * s_vols[2]
    v3pre = Q * (-2 * s_vols[2] + dx)

    vol = v1pre * s_engs[1] + v2pre * s_engs[2] + v3pre * s_engs[3]
    # Quadratic fit
    fit = np.polyfit(s_vols, s_engs, 2)
    lstsq_abc = dumps(fit.tolist())
    bm_lstsq = float(-fit[1]) * ev_a3_to_gpa  # APPROXIMATE

    return bulkmod, vol, bm_lstsq, lstsq_abc
