from typing import Tuple as T
from json import loads, dumps
import numpy as np


def analyze_opt_bulks(
    vols_: str, engs_: str, mags_: str, strains_: str, xcoh: bool,
    xbm: bool, xlat: bool, par: str
) -> T[float, float, str, float, float, bool]:
    '''Analysis on serialized info about 5 jobs centered on optimum.'''
    DEBUG = False
    ev_a3_to_gpa = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    serial = [vols_, engs_, mags_, strains_]
    s_vols, s_engs, s_mags, s_strains = map(loads, serial)
    only_coh = xcoh and not (xbm or xlat)
    if DEBUG:
        if par == 'BaO':
            import pdb
            pdb.set_trace()
    if len(s_vols) == 5:
        dx = s_vols[1]-s_vols[0]
        stencil = (-1*s_engs[4] + 16*s_engs[3] - 30*s_engs[2]
                   + 16*s_engs[1] - s_engs[0])/(12*dx**2)
        mag, energy = s_mags[2], s_engs[2]
        bulkmod = stencil * s_vols[2] * ev_a3_to_gpa  # APPROXIMATE
        # Quadratic fit
        fit = np.polyfit(s_vols, s_engs, 2)
        lstsq_abc = dumps(fit.tolist())
        bm_lstsq = float(-fit[1]) * ev_a3_to_gpa  # APPROXIMATE
        success = True
    elif only_coh and 0 in s_strains:
        zero_ind = s_strains.index(0)
        mag, energy = s_mags[zero_ind], s_engs[zero_ind]
        bulkmod = bm_lstsq = lstsq_abc = None  # type: ignore
        success = True
    else:
        bulkmod = bm_lstsq = lstsq_abc = mag = energy = None  # type: ignore
        success = False

    return bulkmod, bm_lstsq, lstsq_abc, mag, energy, success
