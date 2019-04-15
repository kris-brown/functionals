from typing import List as L, Tuple as T
from json import load,dumps

def parse_fit(pth : str)->T[L[str],L[str],L[int],L[float],L[str],
                            L[str],L[float],L[float],L[float]]:
    from os import walk
    from os.path import join, getmtime, getsize
    from datetime import datetime  # type: ignore
    try:
        files = [join(dirpath, filename)
             for (dirpath, dirs, files) in walk(pth)
             for filename in (dirs + files) if filename =='result.json' and getsize(join(dirpath, filename)) > 0]
        pths,names,pws,ecs,datas,consts,regs,\
            bmws,latws,ts = \
             [],[],[],[],[],[],[],[],[],[]

        for file in files:
            pth = file[:file.rfind('/')]
            with open(join(pth,'metadata.json'),'r') as f: d = load(f)
            p,c = d['params'], d['calc']

            names.append(d['uid']);  pths.append(pth)
            ecs.append(d['calc']['econv']); pws.append(d['calc']['pw'])

            datas.append(dumps([c['fx']]+c['decays']+[c['msb']]))
            consts.append(' '.join(d['params']['c']))
            regs.append(d['params']['reg'])
            bmws.append(d['params']['bmw']); latws.append(d['params']['lw'])
            ts.append(datetime.fromtimestamp(getmtime(file) + 8*3600)) # add PST shift)

        return (names,pths,ts,pws,ecs,datas, # type: ignore
                consts,regs,bmws,latws)
    except:
        import pdb,traceback
        traceback.print_exc();pdb.set_trace(); assert False
