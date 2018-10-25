from typing import Tuple
from json import loads
from dbgen.core.lists import normalize_list,nub
#########################################################################
def countatm(raw:str)->Tuple[int,int,str,str,str,str,str]:
    """Stoichiometry analysis"""
    nonmets  = [1,2,6,7,8,9,10,16,17,18]
    nums     = sorted([a['number'] for a in loads(raw)['atomdata']])
    nnorm    = normalize_list(nums)
    metnums  = [n for n in nums if n not in nonmets]
    uniqnums = nub(nums)
    comp     = {n:nums.count(n)  for n in uniqnums}
    compnorm = {n:nnorm.count(n) for n in uniqnums}
    metcomp  = {n:metnums.count(n) for n in uniqnums if n not in nonmets}
    consts   = str([a['constrained'] for a in loads(raw)['atomdata']])
    return len(nums),len(uniqnums),str(comp),str(compnorm),str(metcomp),\
            str(nums),consts
