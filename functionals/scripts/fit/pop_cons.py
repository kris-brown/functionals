from typing import List as L, Tuple as T


def pop_cons(params_: str) -> T[L[str], L[str], L[str], L[str],
                                L[float], L[str]]:
    '''Prepares constraints for insertion into database'''
    import json
    from functionals.fit.constraint import constlist

    params = [json.loads(x) for x in params_.split('$')]
    abdata = []
    for c in constlist:
        allpoints = set(c.points)
        for p in params:
            allpoints.update(map(tuple, p.get(c.name, [])))  # type: ignore
        abdata.append(str({(s, a): c.fun(s, a).tolist()
                           for s, a in allpoints}))

    name, kind, ps, vs, fs = map(list, zip(*[
        (c.name, c.kind, json.dumps(c.points), c.val,
         c.fun.__name__) for c in constlist]))

    return name, abdata, kind, ps, vs, fs  # type: ignore
