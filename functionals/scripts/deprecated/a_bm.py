from typing import Tuple as T
import numpy as np
from numpy.linalg import inv
from json import dumps, loads


def a_bm(engs_: str,
         vols_: str,
         contribs_: str,
         coefs_: str,
         bm_expt_: float,
         allengs_: str,
         allvols_: str,
         vol: float, xvol: float
         ) -> T[str, float, str, float]:
    '''
    E(Vol) = Ex(Vol) + Enx(Vol) = A * Vol² + B * Vol + C
     - A :: eV/A^-6, B :: eV/A^-3, C :: eV

    # FIRST: BULK MODULUS
    BM_expt = 2 * A * [[eV/A³ to GPa]]

        - [[eV/A³ to GPa]] = (10**-9) * (1.602 * 10**-19) * (10**10)**3
        - GPa = GJ/m³

    Given any e (energy) and vol (volume) vectors we can solve: e = V · x̅
        - x̅ = (Vᵗ·V)⁻¹ · Vᵗ · e
        - x̅ = least square fit coefficients : [[A],[B],[C]]
        - V = a "vandermonde" matrix generated purely by the volume data
            = [[1, vol1, vol1²]]
                    ...
              [[1, vol5, vol5²]]

    so = (Aᵗ·A)⁻¹·Aᵗ·e
    or A = [1,0,0] · (Aᵗ·A)⁻¹·Aᵗ·e

    we just need the energy vector as a function of an arbitrary beef coef vect
    we need to hold the non-x energy constant and add in the variable x part
    E = Enx + C · x̅ = Edft - C·x̅_old + C· x̅

    Making our final expression for BM_expt:
    BM_expt = 2 * CONV * [1,0,0] · [(Vᵗ·V)⁻¹·Vᵗ]· (Edft - C·x̅_old + C· x̅ )
    BM_expt = a_bm · x̅ + b_bm
        - a_bm = 2 * CONV * [1,0,0] · [(Vᵗ·V)⁻¹·Vᵗ]·C
        - b_bm = 2 * CONV * [1,0,0] · [(Vᵗ·V)⁻¹·Vᵗ]·(Edft - C·x̅_old)

    with
        - CONV = (10**-9) * (1.602 * 10**-19) * (10**10)**3
        - V    = "vandermonde" volume matrix
        - Edft = 5 element vector of DFT energies at different strains
        - C    = 5 x 64 matrix with exchange contributions for each strain
        - x̅    = Arbitrary input vector
        - x̅_old= BEEF coefficients used for the fitted calculations

    For volume, note that  minimum of parabola is at - 0.5 * B / A
        - we'll use A = 0.5 * BM_expt
        - our BM_expt is in units GPa, so need to convert to eV/A^3
    V_expt = - ([0,1,0] · [(Vᵗ·V)⁻¹·Vᵗ]·(Edft - C·x̅_old + C· x̅ )) / BM_expt
    V_expt = a_l · x̅ + b_l
     - a_l = - 1/BM_expt * ([0,1,0] · [(Vᵗ·V)⁻¹·Vᵗ] · C
     - b_l = - 1/BM_expt * ([0,1,0] · [(Vᵗ·V)⁻¹·Vᵗ] · (Edft - C·x̅_old)

    '''

    # Common stuff to preprocessing both BM and Lattice data
    # ----------------------------------------------------------
    bm_conv = (10**-9) * (1.602 * 10**-19) * (10**10)**3

    volumes = np.array(loads(vols_))             # 5 element array, A^3
    energies = np.array(loads(engs_))             # 5 element array, eV

    contribs = np.array(loads(contribs_))
    N = len(volumes)  # number of jobs used
    coefs = np.array(loads(coefs_))            # 64 element array

    fit = np.polyfit(loads(allvols_), loads(allengs_), 2)
    if True:  # bm_expt_ is None:
        bm_expt = 2 * fit[0] * volumes[0] * bm_conv  # APPROXIMATE
    else:
        bm_expt = float(bm_expt_)

    c = np.vstack(contribs)  # n points x 64
    e_nonx = energies - c @ coefs  # len-5 vector

    vander = np.vstack((np.ones(N), volumes, volumes**2)).T  # 5 x 3
    vinv = inv(vander.T @ vander)                            # 3 x 3
    solv = vinv @ vander.T                                   # 3 x 5

    bm_base = 2 * bm_conv * volumes[0] * (np.array([0, 0, 1]) @ solv)
    a_bm = bm_base @ c
    b_bm = bm_base @ e_nonx

    bm_expt = bm_base @ c @ coefs + (bm_base @ e_nonx)
    lat_base = -bm_conv * volumes[0] / bm_expt * (np.array([0, 1, 0]) @ solv)
    a_l = lat_base @ c
    b_l = lat_base @ e_nonx

    return dumps(a_bm.tolist()), b_bm, dumps(a_l.tolist()), b_l
