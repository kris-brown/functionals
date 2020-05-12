def parse_procar(pth: str):
    """Parse procar, return (string) error if any"""
    from re import findall, MULTILINE
    import os

    # elem = pth.split('/')[-1]
    ppth = os.path.join(pth, 'PROCAR')
    if not(os.path.exists(ppth)):
        return ppth+" does not exist"
    with open(ppth, 'r') as fi:
        pro = fi.read()

    # from WebElements: note if a whole shell is filled it can be ignored.
    # I.e. Number of S electrons is modulo 2, P is modulo 6, D is modulo 10, etc.
    configs = dict(H=[1,0,0,0], He=[2,0,0,0], Li=[1,0,0,0], Be=[2,0,0,0], B=[2,1,0,0], C=[2,2,0,0], N=[2,3,0,0], O=[2,4,0,0], F=[2,5,0,0], Ne=[2,6,0,0], Na=[1,0,0,0], Mg=[2,0,0,0], Al=[2,1,0,0], Si=[2,2,0,0], P=[2,3,0,0], S=[2,4,0,0], Cl=[2,5,0,0], Ar=[2,6,0,0], K=[1,0,0,0], Ca=[2,0,0,0], Sc=[2,0,1,0], Ti=[2,0,2,0], V=[2,0,3,0], Cr=[1,0,5,0], Mn=[2,0,5,0], Fe=[2,0,6,0], Co=[2,0,7,0], Ni=[2,0,8,0], Cu=[1,0,10,0], Zn=[2,0,10,0], Ga=[2,1,10,0], Ge=[2,2,10,0], As=[2,3,10,0], Se=[2,4,10,0], Br=[2,5,10,0], Kr=[2,6,10,0], Rb=[1,0,0,0], Sr=[2,0,0,0], Y=[2,0,1,0], Zr=[2,0,2,0], Nb=[1,0,4,0], Mo=[1,0,5,0], Tc=[2,0,5,0], Ru=[1,0,7,0], Rh=[1,0,8,0], Pd=[0,0,10,0], Ag=[1,0,10,0], Cd=[2,0,10,0], In=[2,1,10,0], Sn=[2,2,10,0], Sb=[2,3,10,0], Te=[2,4,10,0], I=[2,5,10,0], Xe=[2,6,10,0], Cs=[1,0,0,0], Ba=[2,0,0,0], La=[2,0,1,0], Ce=[2,0,1,1],  Pr=[2,0,0,3],  Nd=[2,0,0,4],  Pm=[2,0,0,5],  Sm=[2,0,0,6],  Eu=[2,0,0,7],  Gd=[2,0,1,7],  Tb=[2,0,0,9],  Dy=[2,0,0,10],  Ho=[2,0,0,11],  Er=[2,0,0,12],  Tm=[2,0,0,13],  Yb=[2,0,0,14],  Lu=[2,0,1,14],  Hf=[2,0,2,14],  Ta=[2,0,3,14],  W=[2,0,4,14],  Re=[2,0,5,14],  Os=[2,0,6,14],  Ir=[2,0,7,14],  Pt=[1,0,9,14],  Au=[1,0,10,14],  Hg=[2,0,10,14],  Tl=[2,1,10,14],  Pb=[2,2,10,14],  Bi=[2,3,10,14],  Po=[2,4,10,14],  At=[2,5,10,14],  Rn=[2,6,10,14])

    # Parse the PROCAR
    pat = r"band\s+\d+\s\#\senergy\s+([-]?\d+.\d+)\s\#\socc.\s+(\d.\d+)\s+ion.*$\s+1\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+"
    vals = [list(map(float, x)) for x in findall(pat, pro, MULTILINE)]
    if not vals:
        pat = r"band\s+\d+\s\#\senergy\s+([-]?\d+.\d+)\s\#\socc.\s+(\d.\d+)\s+ion.*$\s+1\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+"
        vals = [list(map(float, x)) for x in findall(pat, pro, MULTILINE)]
        vals = [v[:5]+[0,v[-1]] for v in vals]

    # Initialize accumulators
    occ_eng, unocc_eng, ss, ps, ds, fs = [], [], [0.], [0.], [0.], [0.]

    # Iterate through parsed PROCAR
    for eng, occ, s, p, d, f, total in vals:
        spdf = [s, p, d, f]
        if occ%1 > 0.02 and occ%1 < 0.98:
            return "all occupations must be 0/1: %f" % occ
        elif occ > 0.98:
            occ_eng.append(eng)
            # ALLOW MIXED STATES
            # rats = [x/total for x in spdf]
            # if sorted(rats)[-2] > 0.1:
            #     return "S/P/D/F "+str(rats)
            # for tot, rat in zip([ss, ps, ds, fs], rats):
            #     tot[0] += rat*occ
        else:
            unocc_eng.append(eng)

    # Make sure LUMO < HOMO
    # if unocc_eng and occ_eng and min(unocc_eng) < max(occ_eng):
    #     return "bad eng ordering: \nunocc %s \nvs occ %s" %(sorted(unocc_eng), sorted(occ_eng))

    # No longer check whether we have correct occupancy
    # full = [2,6,10,14]  # number of electrons in filled shell
    # config = [int(round(xs[0])%l) for l, xs in
    #           zip(full,[ss,ps,ds,fs])]
    # ref = [x%l for l, x in zip(full,configs[elem])]
    # if config != ref:
    #     return "got %s != %s expected"%(config, ref)
    # return ''


if __name__ == '__main__':
    import sys
    print(parse_procar(sys.argv[1]))