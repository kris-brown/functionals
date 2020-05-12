from typing import List as L, Tuple as T, Any


def pareto(Xs: L[Any], Ys: L[Any]) -> T[L[Any], L[Any]]:
    maxX, maxY = False, False
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]
    for pair in myList[1:]:
        if maxY:
            if pair[1] >= p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[1] <= p_front[-1][1]:
                p_front.append(pair)
    p_frontX, p_frontY = zip(*p_front)
    return p_frontX, p_frontY
