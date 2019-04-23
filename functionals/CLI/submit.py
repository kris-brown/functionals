# External modules
from typing      import Any, Set as S, Union as U, List as L, Callable as C
from abc         import ABCMeta, abstractmethod
from os.path     import join,exists
from os          import chdir,system,mkdir,listdir,environ,getenv
from os          import mkdir
from subprocess  import check_output
from random      import random, choice

from ase        import Atoms            # type: ignore
from ase.io     import write,read       # type: ignore
from ase.data   import chemical_symbols # type: ignore

from functionals.CLI.submit_parser import parser
"""
A CLI interface designed to submit atomic and bulk calculations en masse

ex.
python
    functionals/submit.py
    --target=/scratch/users/ksb/atomization/auto/bulks
    --src=/scratch/users/ksb/functionals/data/structures

python
    functionals/submit.py

"""
################################################################################

# Constants
sher    = bool(getenv('SHERLOCK'))
pp_root = '/home/users/ksb/vossj/'  if sher else '/nfs/slac/staas/fs1/g/suncat/ksb/'
ppth    = pp_root + 'vasp/potpaw_PBE'

def readfile(pth:str)->str:
    with open(pth,'r') as f: return f.read()

def get_script(s:str)->str:
    funcroot = environ['FUNCTIONALS_ROOT']
    with open(join(funcroot,'functionals/scripts/'+s),'r') as f:
        return f.read()

def safeMkdir(pth:str) -> None:
    if not exists(pth): mkdir(pth)

# Magnetic moments
emag = {'Ni': 2, 'Rb': 1, 'Pt': 2, 'Ru': 4, 'S': 2, 'Na': 1, 'Nb': 5, 'Mg': 0, 'Li': 1, 'Pb': 2, 'Pd': 0, 'Ti': 2, 'Te': 2, 'Rh': 3, 'Ta': 3, 'Be': 0, 'Ba': 0, 'As': 3, 'Fe': 4, 'Br': 1, 'Sr': 0, 'Mo': 6, 'He': 0, 'C': 2, 'B': 1, 'P': 3, 'F': 1, 'I': 1, 'H': 1, 'K': 1, 'Mn': 5, 'O': 2, 'Ne': 0, 'Kr': 0, 'Si': 2, 'Sn': 2, 'W': 4, 'V': 3, 'Sc': 1, 'N': 3, 'Os': 4, 'Se': 2, 'Zn': 0, 'Co': 3, 'Ag': 1, 'Cl': 1, 'Ca': 0, 'Ir': 3, 'Al': 1, 'Cd': 0, 'Ge': 2, 'Ar': 0, 'Au': 1, 'Zr': 2, 'Ga': 1, 'In': 1, 'Cs': 1, 'Cr': 6, 'Cu': 1, 'Y' : 1, 'Tc' : 5, 'Sb':3,'Xe':0, 'Hf':2, 'Re':5,'Hg':0,'Tl':1}
allelems = [chemical_symbols.index(x) for x in emag.keys()] # all atomic numbers in {emag}
bmag = {'Li_bcc': None, 'Na_bcc': None, 'K_bcc': None, 'Rb_bcc': None, 'Ca_fcc': None, 'Sr_fcc': None, 'Ba_bcc': None, 'V_bcc': None, 'Nb_bcc': None, 'Ta_bcc': None, 'Mo_bcc': None, 'W_bcc': None, 'Fe_bcc': 2.22, 'Rh_fcc': None, 'Ir_fcc': None, 'Ni_fcc': 0.64, 'Pd_fcc': None, 'Pt_fcc': None, 'Cu_fcc': None, 'Ag_fcc': None, 'Au_fcc': None, 'Al_fcc': None, 'C_diamond': None, 'Si_diamond': None, 'Ge_diamond': None, 'Sn_diamond': None, 'Pb_fcc': None, 'Cd_hcp': None, 'Co_hcp': 1.72, 'Os_hcp': None, 'Ru_hcp': None, 'Zn_hcp': None, 'Ti_hcp': None, 'Zr_hcp': None, 'Sc_hcp': None, 'Be_hcp': None, 'Mg_hcp': None, 'LiH_b1': None, 'LiF_b1': None, 'LiCl_b1': None, 'NaF_b1': None, 'NaCl_b1': None, 'MgO_b1': None, 'MgS_b1': None, 'CaO_b1': None, 'TiC_b1': None, 'TiN_b1': None, 'ZrC_b1': None, 'ZrN_b1': None, 'VC_b1': None, 'VN_b1': None, 'NbC_b1': None, 'NbN_b1': None, 'FeAl_b2': 0.35, 'CoAl_b2': None, 'NiAl_b2': None, 'BN_b3': None, 'BP_b3': None, 'BAs_b3': None, 'AlN_b3': None, 'AlP_b3': None, 'AlAs_b3': None, 'GaN_b3': None, 'GaP_b3': None, 'GaAs_b3': None, 'InP_b3': None, 'InAs_b3': None, 'SiC_b3': None, 'KBr_b1': None, 'CaSe_b1': None, 'SeAs_b1': None, 'RbI_b1': None, 'LiI_b1': None, 'CsF_b1': None, 'CsI_b2': None, 'AgF_b1': None, 'AgCl_b1': None, 'AgBr_b1': None, 'CaS_b1': None, 'BaO_b1': None, 'BaSe_b1': None, 'CdO_b1': None, 'MnO_b1': 2.4, 'MnS_b1': 2.4, 'ScC_b1': None, 'CrC_b1': 0.6, 'MnC_b1': 1.2, 'FeC_b1': None, 'CoC_b1': None, 'NiC_b1': None, 'ScN_b1': None, 'CrN_b1': 1.3, 'MnN_b1': 1.6, 'FeN_b1': 1.3, 'CoN_b1': None, 'NiN_b1': None, 'MoC_b1': None, 'RuC_b1': None, 'RhC_b1': None, 'PdC_b1': None, 'MoN_b1': None, 'RuN_b1': None, 'RhN_b1': None, 'PdN_b1': None, 'LaC_b1': None, 'TaC_b1': None, 'WC_b1': None, 'OsC_b1': None, 'IrC_b1': None, 'PtC_b1': None, 'LaN_b1': None, 'TaN_b1': None, 'WN_b1': None, 'OsN_b1': None, 'IrN_b1': None, 'PtN_b1': None}

class Calc(object):
    def __init__(self,
                 xc     : str,
                 sigma  : float = 0.05,
                 econv  : float = 5e-3,
                 magmom : float = None,
                 kpts   : tuple = (1,1,1)
                ) -> None:
        assert xc in ['PBE',"BEEF",'SCAN']

        self.pw    = 900 if xc == 'SCAN' else 1000;
        self.sigma = sigma
        self.econv = econv;
        self.kpts  = kpts
        self.xc    = xc
        self.mag   = magmom

        # Hardcode these for now
        self.a11,self.a12,self.a13,self.a14,self.a15 = 2.,2.5,3.,4.5,6.
        self.msb = 4.

    def magmom(self,n_atoms : int) -> str:
        if not self.mag: return 'ISPIN = 1'
        elif n_atoms==1: return 'ISPIN = 2\nNUPDOWN = %f'%self.mag
        else:            return 'ISPIN = 2\nMAGMOM = '+ ' '.join(n_atoms * [('%s' % self.mag)])

class Job(metaclass = ABCMeta):
    def __init__(self, calc : Calc, time : int = 1) -> None:
        assert time >= 0
        self.calc = calc
        self.time = time

    def __repr__(self)->str:
        return str(self)

    def _submit(self, pth : str) -> None:
        self.singlepoint(pth)
        sunc   = choice(['','','','2','3'])
        cd    = {'':8,'2':12,'3':16} # core-dict
        cores = cd[sunc]
        bash = join(pth,'subVASP.sh')
        sub  = self.subsunc
        with open(bash,'w') as g:
            g.write(sub.format(Job=self))
        chdir(pth);
        system('chmod 755 '+bash)
        if self.calc.xc == 'SCAN':
            cmd = './subVASP.sh'
            args = [] # type: list
        else:
            cmd = 'bsub -n {} -W{}:09 -q suncat{} {}'
            args = [cores,self.time,sunc,bash]
        system(cmd.format(*args))

    def done(self,spth:str)->bool:
        ozpth,outpth = [join(spth,f) for f in ['OSZICAR','OUTCAR']]

        don = exists(ozpth) and (len(readfile(ozpth).split()) < 800)            \
                             and ('General timing' in readfile(outpth))
        return don

    @property
    @abstractmethod
    def subsunc(self)->str: raise NotImplementedError

    @abstractmethod
    def submit(self,
               sher  : bool ,
               retry : bool   = False,
               curr  : S[str] = set()
              ) -> None:
        raise NotImplementedError

    @abstractmethod
    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        raise NotImplementedError

class Atomic(Job):
    """
    Calculation of the energy for an isolated atom
    """
    def __init__(self,
                 elem   : int,
                 xc     : str,
                 time   : int   = 10,
                 sigma  : float = 0.01,
                ) -> None:

        self.elem = elem
        self.sym  = sym = chemical_symbols[elem]
        calc = Calc(xc = xc, sigma=sigma, magmom = emag[self.sym], kpts=(1,1,1))
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Atoms<%d>'%(chemical_symbols[self.elem])

    def submit(self,
               sher  : bool ,
               retry : bool   = False,
               curr  : S[str] = set()
              ) -> None:
        '''Creates VASP input files'''
        root = '/nfs/slac/g/suncatfs/ksb/beefjobs/atoms/%s/'%self.calc.xc.lower()

        pth = join(root,self.sym)
        safeMkdir(pth)
        seen = pth in curr
        done = self.done(pth)
        if retry or (not done and join(pth,'subVASP.sh') not in curr):

            atoms = Atoms(numbers    = [self.elem],
                          cell       = [[10.1,0.2,0.1],
                                         [0.2,10.2,0.01],
                                         [0.5,0.1,10.3]],
                          positions  = [[0.01,0.02,0.03]])

            poscar = join(pth,'POSCAR')
            potcar = join(pth,'POTCAR')
            write(poscar,atoms)
            potcmd = 'cat %s/%s/POTCAR > %s'%(ppth,self.sym,potcar)
            system(potcmd) # write potcar
            self._submit(pth)

    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        if self.calc.xc == 'BEEF':
            files = [('INCAR','PBE_ATOMCAR'),('INCAR2','BEEF_ATOMCAR')]
        elif self.calc.xc=='PBE':
            files = [('INCAR','PBE_CAR')]
        elif self.calc.xc=='SCAN':
            files = [('INCAR','SCAN_CAR')]

        for fname,f in files:
            incar = get_script('vasp/'+f)
            incar_str = incar.format(Calc=self.calc,magmom=self.calc.magmom(n_atoms=1))
            with open(join(pth,fname),'w') as g:
                g.write(incar_str)

        kpts = get_script('vasp/KPOINTS')
        with open(join(pth,'KPOINTS'),'w') as g:
            g.write(kpts.format(*self.calc.kpts))

    @property
    def subsunc(self)->str:
        if self.calc.xc == 'SCAN':
            return get_script('vasp/subVASPscan.sh')
        else:
            return get_script('vasp/subVASPatom.sh')

    @staticmethod
    def enmasse(elems     : L[int], curr : S[str],
                xc        : str   = '',
                time      : int   = 10,
                sigma     : float = 0.01,
                sher      : bool  = True,
                ) -> None:
        """Submit many atomic jobs with the same parameters"""
        es = [Atomic(elem=e,xc=xc,time=time,sigma=sigma) for e in elems]
        for e in es:
            e.submit(sher=sher,curr=curr)

class Bulk(Job):
    """
    Multiple singlepoints at different strains to determine energy and bulkmod
    """
    def __init__(self,
                 bulk   : str,
                 xc     : str,
                 time   : int   = 1,
                 sigma  : float = 0.01,
                 econv  : float = 1e-3,
                 strains: list  = list(range(-5,5)) # -9,10
                ) -> None:
        self.atoms   = read(bulk)
        self.name    = bulk[bulk.rfind('/')+1:bulk.find('.')]
        self.strains = strains
        magmom       = bmag[self.name]
        calc = Calc(xc=xc,sigma = sigma, econv = econv, kpts = (10,10,10), magmom = magmom)
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Bulk<%s>'%self.name

    def submit(self,
               sher  : bool ,
               retry : bool   = False,
               curr  : S[str] = set()
              ) -> None:
        root = '/nfs/slac/g/suncatfs/ksb/beefjobs/bulks/%s/'%self.calc.xc.lower()
        pth = join(root,self.name)
        safeMkdir(pth)

        opt_cell = self.atoms.get_cell()
        for strain in self.strains:
            cell = opt_cell * (1 + strain/100.)
            spth = join(pth,'strain_%d'%strain)
            done = self.done(spth)
            if retry or (not done and join(spth,'subVASP.sh') not in curr):
                safeMkdir(spth)
                strained_atoms = self.atoms.copy()
                strained_atoms.set_cell(cell, scale_atoms = True)
                poscar = join(spth,'POSCAR')
                potcar = join(spth,'POTCAR')
                write(poscar,strained_atoms)
                elems = check_output('head -n 1 '+poscar, encoding='UTF-8',shell=True).split()
                potcar = 'cat ' + ' '.join(['%s/%s/POTCAR'%(ppth,x) for x in elems]) + ' > ' + potcar
                system(potcar) # write potcar

                self._submit(spth)

    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        car   = dict(PBE='PBE_CAR',BEEF='INCAR',SCAN='SCAN_CAR')
        incar = get_script('vasp/' + car[self.calc.xc])
        incar_str = incar.format(Calc=self.calc,magmom=self.calc.magmom(n_atoms = len(self.atoms)))
        with open(join(pth,'INCAR'),'w') as g:
            g.write(incar_str)
        kpts = get_script('vasp/KPOINTS')
        with open(join(pth,'KPOINTS'),'w') as g:
            g.write(kpts.format(*self.calc.kpts))

    @property
    def subsunc(self)->str:
        return get_script('vasp/subVASP%s.sh'%('scan' if self.calc.xc=='SCAN' else ''))

    @staticmethod
    def enmasse(pth    : str, curr : S[str],
                xc     : str,
                sher   : bool  = False,
                time   : int   = 1,
                sigma  : float = 0.01,
                econv  : float = 1e-3,
                lo     : int   = 5,
                hi     : int   = 5,
                retry  : bool  = False
                ) -> None:
        """ Create a bunch of bulks from a directory containing .traj files"""
        bulkpaths = [pth] if pth[-5:]=='.traj' else listdir(pth)
        assert all([p[-5:]=='.traj' for p in bulkpaths])
        bulks = [Bulk(bulk = join(pth,p), time = time, sigma = sigma, econv = econv,
                      xc = xc, strains = list(range(-lo,hi))) for p in bulkpaths]

        for b in bulks:
            b.submit(retry=retry, curr = curr, sher=sher)


#########################################

def main() -> None:
    '''Submit either bulk or element singlepoint calculations'''
    args = parser.parse_args()
    assert bool(args.src) ^ bool(args.elems), "Must be submiting Bulk or Atomic jobs"
    #assert args.target,                       "Need a target location to submit jobs"

    currstr = 'squeue -u ksb -o "%Z"' if sher else 'bjobs -o "job_name" -noheader'

    curr   = set(check_output(currstr, encoding='UTF-8',shell=True).split())
    common = dict(xc = args.xc, curr = curr,
                  time = args.time, sigma = args.sigma )

    if args.elems:
        Atomic.enmasse(**{**common,**dict(sher=sher,elems=args.elems)})
    elif args.src:
        Bulk.enmasse(**{**common,**dict(sher = sher, pth = args.src,
                                        lo = args.lo, hi = args.hi)})
    else:
        raise ValueError("Must be submiting Bulk or Atomic jobs")

if __name__=='__main__':
    main()
