from typing import List as L, Tuple as T
from json import load,dumps
from os.path import join,getsize
from os import listdir

def parse_fit(pth : str)->T[L[str],L[str],L[bool],L[int],L[str],
                            L[str],L[float],L[float],L[float]]:
    pths,names,pws,datas,consts,ces,bms,lcs,done = [],[],[],[],[],[],[],[],[]
    req = set(['x%d.json'%i for i in range(3)])
    for fol in listdir(pth):
        dir = join(pth,fol)
        subfolders = set([x for x in listdir(dir) if getsize(join(dir,x)) > 0])
        with open(join(dir,'metadata.json'),'r') as f: d = load(f)
        p,c = d['params'], d['calc']

        names.append(d['uid']);  pths.append(dir)
        pws.append(int(d['calc']['pw']))

        datas.append(dumps(c['fx']))
        consts.append(' '.join(d['params']['c']))
        ces.append(d['params']['ce_scale']);
        bms.append(d['params']['bm_scale']);
        lcs.append(d['params']['lc_scale']);
        done.append(0 if 'fit.json' not in subfolders else
                   (2 if len(req-subfolders)==0 else 1))
        #ts.append(datetime.fromtimestamp(getmtime(file) + 8*3600)) # add PST shift)

    return (names,pths,done,pws,datas, # type: ignore
            consts,ces,bms,lcs)
