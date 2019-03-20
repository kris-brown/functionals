from typing import List as L, Tuple as T
from json import load,dumps

def parse_fit(pth : str)->T[L[str],L[float],L[str],L[int],L[str],L[int],L[float],
                            L[float],L[str],L[str],L[str],L[str],L[str],L[str],
                            L[str]]:
    from os import walk
    from os.path import join, getmtime
    from datetime import datetime  # type: ignore
    files = [join(dirpath, filename)
             for (dirpath, dirs, files) in walk(pth)
             for filename in (dirs + files) if filename =='result.json']
    pths,names,pws,ecs,datas,decayvals,constdens,consts,regs,dataconsts,\
        bmws,latws,ts,res,rts,a1s,a2s,a3s,a4s,a5s,msbs,steps,decays,nsteps = \
         [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    for file in files:
        pth = file[:file.rfind('/')]
        with open(file,'r') as f:                    r = load(f)
        with open(join(pth,'metadata.json'),'r') as f: d = load(f)
        p,c = d['params'], d['calc']

        names.append(d['uid']); decayvals.append(p['a1']); pths.append(pth)
        decays.append(p['decay']); nsteps.append(len(r[1]))
        ecs.append(d['calc']['econv']); pws.append(d['calc']['pw'])
        datas.append(dumps(d['calc']['fx']))
        constdens.append(d['params']['cd'])
        consts.append(' '.join(d['params']['c']))
        regs.append(d['params']['reg'])
        dataconsts.append(d['params']['dc'])
        bmws.append(d['params']['bmw']); latws.append(d['params']['lw'])
        ts.append(datetime.fromtimestamp(getmtime(file) + 8*3600)) # add PST shift)
        rts.append(r[0]);
        steps.append(dumps(r[1]));
        res.append(dumps(r[1][-1]))
        a1s.append(c['a11']); a2s.append(c['a12'])
        a3s.append(c['a13']); a4s.append(c['a14'])
        a5s.append(c['a15']); msbs.append(c['msb'])

    return (names,decayvals,pths,ts,decays,nsteps,steps,res,rts,pws,ecs,datas,constdens, # type: ignore
            consts,regs,dataconsts,bmws,latws,a1s,a2s,a3s,a4s,a5s,msbs)
