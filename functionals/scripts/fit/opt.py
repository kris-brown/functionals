from json import loads, dumps,load
from typing import List as L, Tuple as T
import numpy as np # type: ignore

def opt(pth:str)->T[str,str]:
    import numpy as np
    with open(pth+'/result.json','r') as fi: steps = load(fi)
    print('\n\n\n\n',pth)
    threshold = 1e-3

    def pareto_frontier(Xs:list, Ys:list, maxX:bool = False, maxY:bool = True)->T[list,list]:
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
        import matplotlib.pyplot as plt # type: ignore
        if len(l[0])!=3:
            import pdb;pdb.set_trace()
        dic = {c:i for i,(_,_,c) in enumerate(l)} # remove duplicates
        _,losses,cviols = zip(*l)
        if len(losses) <= 2: return 0

        cv,lo = pareto_frontier(cviols,losses)
        #import pdb;pdb.set_trace()
        dlo = np.gradient(lo,cv)
        for c,dl in reversed(list(zip(cv,dlo))):
            if dl > threshold:
                break

        #plt.scatter(cv,lo); plt.show()
        #plt.scatter(cv,dlo); plt.show()

        #import pdb;pdb.set_trace()
        return 1 + dic[c]

    inds = [f(x[1:]) for x in steps]
    out = (dumps(inds),dumps([step[ind][0] for step,ind in zip(steps,inds)]))
    return out
