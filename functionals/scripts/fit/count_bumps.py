import numpy as np
from typing import List
from functionals.scripts.fit.legcountzeros import countzeros_in_s_at_alpha
alphas = np.array([0, 0.5, 1, 1.5, 2, 100])


def count_bumps(x: np.ndarray, initial_bumps: List[List[int]]) -> int:

    # how much to penalize bumps in 1st, 2nd, 3rd, etc. derivative
    bmpweight = np.array([10, 5, 3, 1, 0, 0, 0, 0])

    a0 = sum([bmpweight @ np.clip(np.array(bump_at_alpha) - ibump_at_alpha,
                                  0, 100)
              for ibump_at_alpha, bump_at_alpha
              in zip(initial_bumps, all_bumps(x))])
    return a0


def all_bumps(x: np.ndarray) -> List[List[int]]:
    return countzeros_in_s_at_alpha(x, alphas)
