from typing import Tuple as T
import numpy as np  # type: ignore
from json import dumps, loads
from functools import reduce


def a_ce(atom_engs_: str,
         atom_contribs_: str,
         bulk_eng: float,
         bulk_contribs_: str,
         n_atoms: int,
         coefs_: str,
         name: str
         ) -> str:
    '''
    Returns vectors which, when dotted with a BEEF coefficient vector and taking
    into account an offset, will yield a cohesive energy in eV/atom

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
    from functools import reduce
    all_atom_engs = [float(x) for x in atom_engs_.split(',')]
    all_atom_contribs = [np.array(loads(x)) for x in atom_contribs_.split('$')]
    atom_engs = sum(all_atom_engs)
    atom_contribs = reduce(np.add, all_atom_contribs)

    # just take the min energy structure's contribs
    bulk_contribs = np.array(loads(bulk_contribs_)[0])

    coefs = np.array(loads(coefs_))  # calculator coefs used in the DFT
    dC = atom_contribs - bulk_contribs
    dE = atom_engs - float(bulk_eng)
    eNonxc = (dE - (dC @ coefs)) / n_atoms

    X = dC / n_atoms

    ms2 = np.array([1.2909350231501469, 0.2465804012284167, -0.03821209014569475, 0.00529705394712429, -0.0007238997245930972, 5.328792422471612e-05, -5.1621942345037315e-05, -4.815694813585898e-05, 0.03151729199875705, -0.05207745850199544, 0.024414808963883195, -0.004361104048983925, 0.0006958550589244329, 1.23444049029741e-05, 0.00011675664799358277, 0.00012931051714575557, -3.984507626043427e-05, 2.0741293699137517e-05, 4.910807626956804e-05, -0.00011955322579136519, -5.192392243427056e-05, -0.00012842806378204823, -0.00013417985726866814, -0.00016678490487096736, 3.6804163267903376e-05, -1.917776957893514e-05, -4.540134087044211e-05, 0.00011033504367016455, 4.820973545648982e-05, 0.00011878859525861501, 0.00012397649516404196, 0.00015431369857345442, -
                    2.6467618673671896e-05, 1.3812919294735846e-05, 3.269598858826801e-05, -7.924699408094865e-05, -3.493487062077754e-05, -8.559944901719986e-05, -8.920008924544498e-05, -0.00011125176685256548, 1.4614502947494727e-05, -7.640968977915107e-06, -1.8083052316753068e-05, 4.3690618400164006e-05, 1.9466090227040113e-05, 4.738006039082977e-05, 4.928000817903369e-05, 6.161257946776035e-05, -6.048969514884463e-06, 3.1746395247115385e-06, 7.510695728058644e-06, -1.8034616473329203e-05, -8.18734544611087e-06, -1.9696146222487697e-05, -2.0423598530726044e-05, -2.5643594612769123e-05, 1.544954038475419e-06, -8.163516253199092e-07, -1.9309474669922123e-06, 4.586015422460599e-06, 2.145156585619744e-06, 5.066002704462669e-06, 5.230302978950324e-06, 6.6110885027564675e-06])

    return dumps([X.tolist(), float(eNonxc)])
