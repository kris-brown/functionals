from typing import List as L, Tuple as T, Optional as O

def parse_const(pth:str) -> T[L[O[str]],L[O[str]],L[O[float]],L[O[float]],L[O[float]],L[O[str]]]:
    from os.path import join
    from os      import environ
    from csv     import reader
    n_,d_,s_,a_,v_,t_ = [],[],[],[],[],[]
    with open(pth,'r') as f:
        r = reader(f); next(r)
        for name,desc,s,a,v,t in r:
            s_.append(float(s) if s else None);
            a_.append(float(a) if a else None);
            v_.append(float(v) if v else None);
            n_.append(name or None); d_.append(desc or None); t_.append(t or None)
    return n_,d_,s_,a_,v_,t_
