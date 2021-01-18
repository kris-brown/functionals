from typing import Tuple as T
import numpy as np # type: ignore
from json import dumps, loads
from functools import reduce

def a_ce(atom_engs_     : str,
         atom_contribs_ : str,
         bulk_engs_     : str,
         bulk_contribs_ : str,
         n_atoms        : int,
         coefs_          : str,
         ) -> T[str,str]:
    '''
    Returns vectors which, when dotted with a BEEF coefficient vector and taking
    into account an offset, will yield a cohesive energy in eV
    (one vector+offset pair for each value of a1)

    CE_expt = Ea - Eb
            = ((Exa+Enxa) - (Exb+Enxb)) / n_atoms
            = (∆Ex + ∆Enx) / n_atoms

    - CE_expt : Experimental cohesive energy
    - x       : exchange
    - nx      : non-exchange
    - a       : atomic
    - b       : bulk
    - ∆???    : ???a - ???b

    Furthermore we have:
        Ex   = Contribs ∙ x̅
        Edft = Ex + Enx
        Enx  = Edft - C ∙ old_x̅

    ∆Ex = ∆C ∙ x̅
    ∆Enx= (Edfta - (Ca ∙ old_x̅)) - (Edftb - (Cb ∙ old_x̅))
        = Edfta - Edftb + (Cb ∙ old_x̅) - (Ca ∙ old_x̅)
        = ∆Edfta + (∆C ∙ old_x̅)

    Thus, for the optimal x̅ coefficient array, we have:
        CE_expt = ∆C/n_atoms ∙ x̅ + ∆Enx/n_atoms
    '''

    try:

        all_atom_engs     = [float(x) for x in atom_engs_.split(',')]
        all_atom_contribs = [np.array(loads(x)) for x in atom_contribs_.split('$')]
        atom_engs     = sum(all_atom_engs)
        atom_contribs = reduce(np.add,all_atom_contribs)
        bulk_contribs = np.array(loads(bulk_contribs_)[0]) # just take the min energy structure's contribs
        bulk_eng      = float(loads(bulk_engs_)[0]) # just take the min energy
        coefs         = np.array(loads(coefs_))
    except Exception as e:
        print(e)
        ;breakpoint()
        assert False
    dC = atom_contribs - bulk_contribs
    dE = atom_engs - bulk_eng
    eNonxc = (dE - (dC @ coefs)) / n_atoms

    X = dC / n_atoms
    return dumps([x.tolist() for x in X]),dumps(eNonxc.tolist())
