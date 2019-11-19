from typing import List as L, Tuple as T
from json import load,dumps

def parse_atoms(root:str)->T[L[str],L[str],L[int],L[int],L[str],L[str],L[float],L[bool],L[int],L[int],L[float]]:
    from os      import listdir, environ
    from os.path import join, exists
    from re      import findall,compile,MULTILINE
    from functionals.scripts.load.parse_incar import parse_incar
    from functionals.scripts.load.parse_contribs_vasp import parse_contribs_vasp
    from functionals.scripts.load.parse_eigenval import parse_eigenval
    bf = '/Users/ksb/functionals/data/beefs/beef.json'
    symbols = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
               'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
               'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
               'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
               'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
               'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
               'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
               'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
               'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
               'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
               'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
               'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    pths,names,pws,fxs,contribs,engs,intoccs,nums,tmags,mags,rts \
        = [],[],[],[],[],[],[],[],[],[],[]
    # try:
    # Traverse directory hierarchy with 3 levels (overall,batch,element)
    stordirs = listdir(root)
    for sd in stordirs: # 'date_description'
        assert not sd[:3] == 'old', sd
        folder = join(root,sd)
        for atom in listdir(folder):
            pth = join(folder,atom)

            # Verify job completed successfully
            if exists(join(pth,'OUTCAR')):
                with open(join(pth,'OUTCAR'),'r')  as f: outcar = f.read()
                with open(join(pth,'OSZICAR'),'r') as f: nsteps = sum(1 for _ in f)
            else: outcar =  ''

            if 'General timing' in outcar and nsteps < 800: # completed

                # GET INCAR PARAMETERS
                incar = parse_incar(join(pth,'INCAR'))
                pws.append(int(incar['encut']));
                tmags.append(int(incar.get('magmom') or 0))

                # GET XC-SPECIFIC INFO
                if incar['metagga'] == 'BF':
                    a1msb = [incar[x] for x in ['a1'+y for y in '12345']+['msb']]
                    with open(bf, 'r') as f: beef = load(f)
                    fxs.append(dumps([beef]+a1msb))
                    contribs.append(parse_contribs_vasp(outcar))
                else:
                    contribs.append(None) # type: ignore
                    if incar['metagga'] == 'SCAN':
                        fxs.append("SCAN")
                    elif incar['metagga'] is None and incar['gga']=='PE':
                        fxs.append('PBE')
                    else:
                        import pdb;pdb.set_trace()
                        raise ValueError()

                # GET ENERGY
                pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
                match = findall(pat, outcar); assert match
                engs.append(float(match[-1]))
                # Get runtime
                patt = r'Elapsed time \(sec\):\s+(\d+\.\d+)'
                match = findall(patt, outcar); assert match
                rts.append(int(float(match[-1])/60.))

                # GET OCCUPATION NUMBERS
                with open(join(pth,'EIGENVAL'),'r') as g:
                    io,mag = parse_eigenval(g.read())
                    intoccs.append(io); mags.append(mag)

                # NAME AND PATH
                num = symbols.index(atom)
                pths.append(pth); names.append(atom); nums.append(num)
    return pths,names,rts,pws,fxs,contribs,engs,intoccs,nums,tmags,mags
    # except Exception as e:
    #     import traceback,pdb
    #     traceback.print_exc(); print(pth,e);pdb.set_trace(); assert False
