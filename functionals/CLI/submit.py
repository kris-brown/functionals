from typing import Set, List, Dict, Any, Optional, Iterable
import argparse
import glob
import math
import os
import os.path
import subprocess
import ase # type: ignore
import ase.data # type: ignore
import ase.io # type: ignore
import functionals.templates   as temp

########
# DATA #
########
ppth = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/potpaw_PBE'
root = '/nfs/slac/g/suncatfs/ksb/beefjobs/'
exe_root = '/nfs/slac/staas/fs1/g/suncat/ksb/vasp/olsherlock/vasp.5.4.1beefcar/bin/vasp_'
struct_root = '/'+os.path.join(*__file__.split('/')[:-3],'data/structures')

# Atomic data
n_electrons = dict(Ag=11, Al=3, As=5, Au=11, B=3, Ba=10, Be=2, Br=7, C=4,
     Ca=10, Cd=12, Cl=7, Co=9, Cr=6, Cs=9, Cu=11, F=7, Fe=8, Ga=3,
     Ge=4, H=1, He=2, Hf=4, Hg=12, I=7, In=3, Ir=9, K=7, Kr=8, La=11,
     Li=1, Mg=2, Mn=7, Mo=6, N=5, Na=1, Nb=13, Ne=8, Ni=10, O=6,
     Os=8, P=5, Pb=4, Pd=10, Pt=10, Rb=9, Re=7, Rh=9, Ru=8, S=6,
     Sb=5, Sc=3, Se=6, Si=4, Sn=4, Sr=10, Ta=5, Te=6, Ti=4,
     Tl=3, V=5, W=6, Xe=8, Y=11, Zn=12, Zr=12)
atommag = dict(Ni= 2, Rb= 1, Pt= 2, Ru= 4, S= 2, Na= 1, Nb= 5, Mg= 0,
        Li= 1, Pb= 2, Pd= 0, Ti= 2, Te= 2, Rh= 3, Ta= 3, Be= 0,
        Ba= 0, As= 3, Fe= 4, Br= 1, Sr= 0, Mo= 6, He= 0, C= 2,
        B= 1, P= 3, F= 1, I= 1, H= 1, K= 1, Mn= 5, O= 2, Ne= 0,
        Kr= 0, Si= 2, Sn= 2, W= 4, V= 3, Sc= 1, N= 3, Os= 4,
        Se= 2, Zn= 0, Co= 3, Ag= 1, Cl= 1, Ca= 0, Ir= 3, Al= 1,
        Cd= 0, Ge= 2, Au= 1, Zr= 2, Ga= 1, In= 1, Cs= 1,
        Cr= 6, Cu= 1, Y = 1, Sb=3,Xe=0, Hf=2, Re=5,
        Hg=0,Tl=1,La=0)

# Bulk data
bulkmag = dict(Li= None, Na= None, K= None, Rb= None,
    Ca= None, Sr= None, Ba= None, V= None,
    Nb= None, Ta= None, Mo= None, W= None,
    Fe= 2.22, Rh= None, Ir= None, Ni= 0.64,
    Pd= None, Pt= None, Cu= None, Ag= None,
    Au= None, Al= None, C= None, Si= None,
    Ge= None, Sn= None, Pb= None, Cd= None,
    Co= 1.72, Os= None, Ru= None, Zn= None, Ti= None,
    Zr= None, Sc= None, Be= None, Mg= None, LiH= None,
    LiF= None, LiCl= None, NaF= None, NaCl= None, MgO= None,
    MgS= None, CaO= None, TiC= None, TiN= None, ZrC= None,
    ZrN= None, VC= None, VN= None, NbC= None, NbN= None,
    FeAl= 0.35, CoAl= None, NiAl= None, BN= None, BP= None,
    BAs= None, AlN= None, AlP= None, AlAs= None, GaN= None,
    GaP= None, GaAs= None, InP= None, InAs= None, SiC= None,
    KBr= None, CaSe= None, SeAs= None, RbI= None, LiI= None,
    CsF= None, CsI= None, AgF= None, AgCl= None, AgBr= None,
    CaS= None, BaO= None, BaSe= None, CdO= None, MnO= 2.4,
    MnS= 2.4, ScC= None, CrC= 0.6, MnC= 1.2, FeC= None,
    CoC= None, NiC= None, ScN= None, CrN= 1.3, MnN= 1.6,
    FeN= 1.3, CoN= None, NiN= None, MoC= None, RuC= None,
    RhC= None, PdC= None, MoN= None, RuN= None, RhN= None,
    PdN= None, LaC= None, TaC= None, WC= None, OsC= None,
    IrC= None, PtC= None, LaN= None, TaN= None, WN= None,
    OsN= None, IrN= None, PtN= None)

appx_bulkmod = dict( # FROM EXPT
    C=454.7,Si=101.3,Ge=79.4,Sn=42.8,SiC=229.1,BN=410.2,BP=168,
    AlN=206,AlP=87.4,AlAs=75,GaN=213.7,GaP=89.6,GaAs=76.7,InP=72,InAs=58.6,
    InSb=46.1,LiH=40.1,LiF=76.3,LiCl=38.7,NaF=53.1,NaCl=27.6,MgO=169.8,Li=13.1,
    Na=7.9,Al=77.1,K=3.8,Ca=15.9,Rb=3.6,Sr=12,Cs=2.3,Ba=10.6,V=165.8,Ni=192.5,
    Cu=144.3,Nb=173.2,Mo=276.2,Rh=277.1,Pd=187.2,Ag=105.7,Ta=202.7,W=327.5,
    Ir=362.2,Pt=285.5,Au=182,ZrC=207.0,Ru=320.8,Be=100.3,Ti=105.1,FeAl=136.1,
    Cd=46.7,Fe=168.3,CaO=115.0,Pb=48.8,Zr=83.3,Mg=35.4,MgS=80.0,NbN=287.0,
    VC=303.0,NbC=300.0,Zn=59.8,TiN=277.0,Co=191.4,CoAl=162.0,BAs=138.0,TiC=233.0,
    NiAl=166.0,VN=268.0,Os=418.0,Sc=43.5,ZrN=240.0,
    # FROM SCAN
    AgBr=47.0,AgF=77.4,BaO=74.3,CdO=161.4,
    CoC=340.3,FeC=359.3,FeN=194.4,LaC=83.9,
    IrN=338.0,IrC=351.7,KBr=17.6,LaN=128.4,
    MnC=256.9,MnS=73.8,MoC=367.1,NiC=292.8,
    OsN=353.1,OsC=377.7,PdN=234.8,RhN=323.7,
    RuC=339.1,ScN=248.8,SeAs=80.5,TaN=376.8,
    WN=401.7,LiI=22.3,BaSe=38.6,MoN=363.9,
    MnN=171.1,PtN=273.8,PdC=244.1,NiN=296.6,
    CsF=25.9,ScC=173.7,CrN=272.6,WC=402.6,
    TaC=355.0,RhC=313.8,CoN=344.4,PtC=297.7,
    MnO=165.4,AgCl=56.1,CrC=294.9,CaSe=56.1,
    CsI=12.896392,CaS=62.919975,
    # PBE/Internet
    RbI=11,RuN=304)

xcdict = dict(pbe = dict(gga='PE'), scan=dict(metagga='SCAN'),
              beef=dict(metagga='BF',lbeefens=True,a11=4.9479, a12=4.9479,
                        a13=4.9479, a14=4.9479, a15=4.9479, msb=1.)) # type: Dict[str,Dict[str,Any]]
################################################################################
# HELPER FUNCTIONS #
def str2bool(v:str)->bool:  return v.lower() in ("yes", "true", "t", "1")

def onezero(one:int ,zero:int) -> str:
    return ('%d*1.0'%one if one else '') + ' %d*0.0'%zero # Render FERWE/FERDO

def get_strain(mat:str, strain_int: int) -> float:
    '''Gives the *lattice* strain, in order to get a *volume* strain.'''
    spacing = 0.25/(appx_bulkmod[mat])**(1./3.)
    return (1 + strain_int * spacing)**(1/3)

def current_jobs() -> Set[str]:
    cmd = 'bjobs -o "job_name" -noheader'
    return set(subprocess.check_output(cmd, encoding='UTF-8',shell=True).split())

def readfile(pth:str)->str:
    """Convert filepath to string of contents."""
    with open(pth,'r') as f: return f.read()

def done(spth : str) -> bool:
    '''Identify if a job has already successfully completed.'''
    ozpth,outpth = [os.path.join(spth,f) for f in ['OSZICAR','OUTCAR']]

    don = os.path.exists(ozpth) and (len(readfile(ozpth).split('\n')) < 800)            \
                         and ('General timing' in readfile(outpth))
    return don

def mk_incar(d: dict, pth: str) -> None:
    """Convert dictionary to INCAR file."""
    def fmt(val: Any) -> str:
        if isinstance(val, bool): return '.TRUE.' if val else '.FALSE.'
        else: return str(val)

    with open(pth, 'w') as f:
        lines = ['{} = {}'.format(k.upper(), fmt(v)) for k, v in d.items()]
        f.write('\n'.join(lines)+'\n')
################################################################################
# MAIN FUNCTIONS #
def submit_atom(name:str, time: int, xc: str, retry: bool, incar: Dict[str,Any],
                curr:Set[str], orbitals: Optional[bool], sunc: str) -> None:
    elem = ase.data.chemical_symbols.index(name)
    atoms = ase.Atoms(numbers=[elem], positions=[[1.,2.,3.]],
                      cell=[[19.,2.,1.], [2.,20,1.], [5.,1,21.]])
    magmom, elec = atommag[name], n_electrons[name]
    nb = max(8, int(math.ceil(0.6*elec+magmom if magmom else elec/2+0.5*1)))
    ferwedo = orbitals if (orbitals is not None) else (elem > 56)
    if ferwedo and magmom:
        nup, ndn = (elec + magmom)//2, (elec - magmom)//2
        magdic = dict(ismear = -2, ferwe = onezero(nup, nb - nup),
                      ferdo = onezero(ndn, nb - ndn))
        #print(orbitals, elem, magmom); import pdb;pdb.set_trace()
    else:
        magdic = dict(sigma=0.0001, nupdown=magmom, ispin=2 if magmom else 1)

    magcar = dict(
        ldiag = False, isym = -1, amix=0.1 if xc=='scan' else 0.2,
        bmix = 0.0001, nbands = nb, **magdic)

    # Submit
    pth = os.path.join(root, 'atoms', xc, name)
    os.makedirs(pth, exist_ok=True)
    notdone = not done(pth)
    notcurr = os.path.join(pth,'subVASP.sh') not in curr
    if retry or (notdone and notcurr):
        _submit(pth,time,atoms,{**incar,**magcar},sunc)

def submit_bulk(mat:str, time: int, xc: str, strains: Iterable[int],
                retry: bool, incar: Dict[str,Any], curr:Set[str], sunc: str,
                ) -> None:
    # Get unstrained traj
    if mat+'.cif' in os.listdir(struct_root):
        orig_atoms = ase.io.read(os.path.join(struct_root, mat+('.cif')))
    else:
        traj = glob.glob(os.path.join(struct_root,mat+('_*.traj')))[0]
        orig_atoms = ase.io.read(traj)
    # Bulk-specific INCAR settings
    magmom = bulkmag[mat]
    magdic = {'ispin':2, 'magmom': ' '.join([str(magmom)]*len(orig_atoms))} \
             if magmom else dict(ispin=1)
    magcar = dict(sigma=0.01, ismear=0, **magdic)

    for strain_int in strains:
        # Strained Traj
        lat_strain = get_strain(mat,strain_int)
        atoms = orig_atoms.copy()
        atoms.set_cell(atoms.get_cell()*lat_strain,scale_atoms=True)
        # Determine whether to submit again or not
        pth = os.path.join(root, 'bulks', xc, mat, 'strain_%d'%strain_int)
        os.makedirs(pth, exist_ok=True) 
        notdone = not done(pth)
        notcurr = os.path.join(pth,'subVASP.sh') not in curr
        if retry or (notdone and notcurr):
            _submit(pth, time, atoms, {**incar,**magcar}, sunc)

def _submit(pth: str, time: int, atoms:ase.Atoms, incar: Dict[str,Any],
            sunc: str) -> None:
    '''General submission for atomic or bulk (at a strain).'''
    # Copy PBE wavecar if possible
    if 'pbe' not in pth:
        pbejob = pth.replace('/beef/','/pbe/').replace('/scan/','/pbe/')
        if done(pbejob):
            os.system('cp {}/WAVECAR {}/WAVECAR'.format(pbejob,pth))

    # Write INCAR
    mk_incar(incar, os.path.join(pth,'INCAR'))

    # Determine executable and kpoints
    scan = incar.get('metagga')=='SCAN'
    if 'atoms' in pth:
        exe, kpts = exe_root + 'gam', (1,1,1)
    else:
        exe, kpts = exe_root + 'std', (10,10,10)
    if scan: exe = 'vasp'

    # KPOINTS
    with open(os.path.join(pth, 'KPOINTS'), 'w') as g:
        g.write(temp.jinja_env.get_template('KPOINTS.jinja').render(kpts=kpts))

    # Write POSCAR and POTCAR
    poscar, potcar = [os.path.join(pth, x) for x in ['POSCAR', 'POTCAR']]
    ase.io.write(poscar,atoms)
    elems = subprocess.check_output('head -n 1 '+poscar, encoding='UTF-8',shell=True).split()
    potcar = 'cat ' + ' '.join(['%s/%s/POTCAR'%(ppth,x) for x in elems]) + ' > ' + potcar
    os.system(potcar)

    # Submission script
    bash = os.path.join(pth,'subVASP.sh')
    with open(bash,'w') as g:
        g.write(temp.jinja_env.get_template('subVASP.jinja').render(exe=exe,scan=scan))
    os.chdir(pth); os.system('chmod 755 '+bash)
    args = [16 if sunc else 8, time, sunc, bash]
    cmd = 'bsub -n {} -W{}:09 -q suncat{} {}'.format(*args)
    # print(cmd);import pdb;pdb.set_trace()
    os.system(cmd)

################################################################################
parser = argparse.ArgumentParser(description  = 'Submit some jobs', allow_abbrev = True)
parser.add_argument('--time', default=10, type=int, help='Walltime')
parser.add_argument('--orb', default=None, type=str2bool, help='Use explicit orbital occupations')
parser.add_argument('--retry', default=False, type=bool, help='Redo jobs')
parser.add_argument('--sunc', default='', type=str, help='Redo jobs')
parser.add_argument('--strain', default = None, help='Two integers for range()',
                    type = lambda x:[int(y) for y in x.split()])
parser.add_argument('--mat', default = '', help='Material name',
                    type = lambda x: list(bulkmag.keys()) if x=='all' else x.split())
parser.add_argument('--elems', default='', help='"all" or list of positive integers',
                    type = lambda x: list(atommag.keys()) if x=='all' else x.split())
parser.add_argument('--xc', help='Copies a file into the working directory as BEEFoftheDay.txt')

def main() -> None:
    '''Submit either bulk or element singlepoint calculations.'''
    # Get/validate command line args
    args = parser.parse_args()
    assert args.xc in ['pbe','scan','beef']
    assert args.sunc in ['','3']
    assert args.time > 0
    curr = current_jobs()

    incar = dict(encut=900., ediff=1e-5, algo='A', istart=2,
                 addgrid=True, ibrion=-1, npar=1, nelm=800, lasph=True,
                 lcharg=True, lwave=True, prec='ACCURATE', **xcdict[args.xc])

    common = dict(incar=incar, time=args.time, xc=args.xc, curr=curr,
                  retry=args.retry, sunc=args.sunc)
    if args.mat:
        assert not args.elems
        for m in args.mat:
            strains = range(*(args.strain or [-4,6]))
            submit_bulk(mat=m, strains=strains, **common)
    else:
        assert args.elems
        for e in args.elems:
            submit_atom(name=e, orbitals=args.orb, **common)


if __name__=='__main__':
    main()
