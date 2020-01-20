from typing import Any, List as L, Tuple as T
import numpy as np


def pop_fitp() -> T[L[float], L[float], L[float], L[str]]:
    import itertools
    import json
    # import functionals.fit.constraint as constr

    def lp(*args: Any) -> L[Any]:
        return list(itertools.product(*args))

    def prod(xl: float, xh: float, nx: int, yl: float, yh: float, ny: int
             ) -> L[Any]:
        return lp(np.linspace(xl, xh, nx),  np.linspace(yl, yh, ny))
    # lowa = prod(0, 4, 5, 0, 0.75, 10,)
    # strict = prod(0.1, 4., 100, 0.1, 4., 30)
    curvs = ['acurvpos', 'acurvneg', 'curvpos', 'curvneg']
    consts = [  # {c.name: [] for c in constr.constlist[3:]},
        dict(),
        dict(**{x: [] for x in curvs}),
        dict(**{x: [] for x in curvs+['zerocurv_hi',
                                      'zerocurv_lo', 'lda_hi', 'lda_lo']})
        #   dict(acurvpos=lowa, acurvneg=lowa),
        #   dict(acurvpos=strict, acurvneg=strict, curvpos=strict,
        #        curvneg=strict),
    ]

    spec = dict(
        consts=list(map(json.dumps, consts)),
        ce_scale=[1000],
        bm_scale=[3.],
        lc_scale=[0.4, 4])
    params = list(itertools.product(*[
        set(x) for x in spec.values()]))  # type: ignore
    con, ce, bm, lc = map(list, zip(*params))
    return ce, bm, lc, con  # type: ignore


if __name__ == '__main__':
    for x in pop_fitp():
        print(x)
