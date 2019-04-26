from json import loads, dumps,load
from typing import List as L, Tuple as T
import numpy as np # type: ignore

def opt(pth:str)->T[str,str]:
    np.warnings.filterwarnings('ignore')

    with open(pth+'/result.json','r') as fi: steps = load(fi)
    threshold = 1e-1
    startup   = 10
    def pareto_frontier(Xs:list, Ys:list, maxX:bool = False, maxY:bool = True) -> T[list,list]:
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse = maxX)
        p_front = [myList[0]]
        for pair in myList[1:]:
            if maxY:
                if pair[1] >= p_front[-1][1]: p_front.append(pair)
            else:
                if pair[1] <= p_front[-1][1]: p_front.append(pair)
        p_frontX,p_frontY = zip(*p_front)
        return p_frontX, p_frontY

    def f(l : L[T[list,float,float]])->int:
        if len(l) <= startup: return 0
        #import matplotlib.pyplot as plt # type: ignore
        dic = {c:startup+i for i,(_,_,c) in enumerate(l[startup:])} # remove duplicates
        _,losses,cviols = zip(*l[startup:])

        cv,lo = pareto_frontier(cviols,losses)
        if len(cv)==1: return dic[cv[0]]

        dlo   = np.gradient(lo,cv)

        for c,dl in reversed(list(zip(cv,dlo))):
            if dl > threshold: break

        #plt.scatter(cv,lo); plt.show();plt.scatter(cv,dlo); plt.show();import pdb;pdb.set_trace()
        return dic[c]

    inds = [f(x[1:]) for x in steps]
    out  = (dumps(inds),dumps([step[ind][0] for step,ind in zip(steps,inds)]))
    return out
