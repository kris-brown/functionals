from typing import Any, List as L, Tuple as T


def pop_fitp() -> T[L[float], L[float], L[float]]:
    import itertools
    # import json
    # import functionals.fit.constraint as constr

    def lp(*args: Any) -> L[Any]:
        return list(itertools.product(*args))

    spec = dict(
        # consts=[dict(pos=3, liebox=3)],
        ce_scale=[1, 4, 10],
        bm_scale=[1, 3],
        lc_scale=[1, ])
    params = list(itertools.product(*spec.values()))
    ce, bm, lc = map(list, zip(*params))
    # con = [json.dumps(c) for c in con_]
    # ab_ = map(constr.Constraint.AB, con_)  # type: ignore
    # ab = [json.dumps([x.tolist(), y.tolist()]) for x, y in ab_]
    return ce, bm, lc  # type: ignore


if __name__ == '__main__':
    for x in pop_fitp():
        print(x)
