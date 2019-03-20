# External modules
from typing      import Any, Set as S, Union as U, List as L, Callable as C
from abc         import ABCMeta, abstractmethod
from os.path     import join,exists
from os          import chdir,system,mkdir,listdir,environ,getenv
from os          import mkdir
from subprocess  import check_output
from argparse    import ArgumentParser
from random      import random, choice

from ase        import Atoms            # type: ignore
from ase.io     import write,read       # type: ignore
from ase.data   import chemical_symbols # type: ignore

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

# Atomic magnetic moments
emag     = {'Ni': 2, 'Rb': 1, 'Pt': 2, 'Ru': 4, 'S': 2, 'Na': 1, 'Nb': 5, 'Mg': 0, 'Li': 1, 'Pb': 2, 'Pd': 0, 'Ti': 2, 'Te': 2, 'Rh': 3, 'Ta': 3, 'Be': 0, 'Ba': 0, 'As': 3, 'Fe': 4, 'Br': 1, 'Sr': 0, 'Mo': 6, 'He': 0, 'C': 2, 'B': 1, 'P': 3, 'F': 1, 'I': 1, 'H': 1, 'K': 1, 'Mn': 5, 'O': 2, 'Ne': 0, 'Kr': 0, 'Si': 2, 'Sn': 2, 'W': 4, 'V': 3, 'Sc': 1, 'N': 3, 'Os': 4, 'Se': 2, 'Zn': 0, 'Co': 3, 'Ag': 1, 'Cl': 1, 'Ca': 0, 'Ir': 3, 'Al': 1, 'Cd': 0, 'Ge': 2, 'Ar': 0, 'Au': 1, 'Zr': 2, 'Ga': 1, 'In': 1, 'Cs': 1, 'Cr': 6, 'Cu': 1, 'Y' : 1, 'Tc' : 5, 'Sb':3,'Xe':0, 'Hf':2, 'Re':5,'Hg':0,'Tl':1}

allelems = [chemical_symbols.index(x) for x in emag.keys()] # all atomic numbers in {emag}

class Calc(object):
    def __init__(self,
                 pw     : int   = 1000,
                 sigma  : float = 0.05,
                 econv  : float = 5e-3,
                 magmom : int   = None,
                 kpts   : tuple = (1,1,1)
                ) -> None:

        self.pw    = pw;
        self.sigma = sigma
        self.econv = econv;
        self.magmom= magmom
        self.kpts  = kpts
        # Hardcode these for now
        self.a11   = 2.
        self.a12   = 2.5
        self.a13   = 3.
        self.a14   = 4.5
        self.a15   = 6.
        self.msb   = 4.


class Job(metaclass = ABCMeta):
    def __init__(self, calc : Calc, time : int = 1) -> None:
        assert time >= 0
        self.calc = calc
        self.time = time

    def __repr__(self)->str:
        return str(self)

    def _submit(self, pth : str) -> None:
        self.singlepoint(pth)
        old = random() < 0.25
        for fi in ['subVASP','subVASP_old']:
            bash = join(pth,fi+'.sh')
            sub  = get_script('vasp/%s.sh'%fi)
            with open(bash,'w') as g:
                g.write(sub.format(Job=self))

            if old == ('old' in fi):
                chdir(pth);
                system('chmod 755 '+bash)
                system('sbatch '+bash)

    def _submit_sunc(self, pth : str) -> None:
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
        cmd = 'bsub -n {} -W{}:09 -q suncat{} {}'
        args = [cores,self.time,sunc,bash]
        x = cmd.format(*args)
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
               pth_  : str,
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
                 time   : int   = 10,
                 sigma  : float = 0.01,
                ) -> None:

        self.elem = elem
        self.sym  = sym = chemical_symbols[elem]
        calc = Calc(sigma=sigma, magmom = emag[self.sym], kpts=(1,1,1))
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Atoms<%d>'%(chemical_symbols[self.elem])

    def submit(self,
               pth_  : str,
               sher  : bool ,
               retry : bool   = False,
               curr  : S[str] = set()
              ) -> None:
        '''Creates VASP input files'''

        pth = join(pth_,self.sym)
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

            if sher:
                self._submit(pth)
            else:
                self._submit_sunc(pth)

    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        for fname,f in [('INCAR','vasp/PBE_ATOMCAR'),('INCAR2','vasp/BEEF_ATOMCAR')]:
            incar = get_script(f)
            incar_str = incar.format(Calc=self.calc)
            if self.calc.magmom is not None:
                incar_str+='\nMAGMOM = '+str(self.calc.magmom)
            with open(join(pth,fname),'w') as g:
                g.write(incar_str)

        kpts = get_script('vasp/KPOINTS')
        with open(join(pth,'KPOINTS'),'w') as g:
            g.write(kpts.format(*self.calc.kpts))

    @property
    def subsunc(self)->str:
        return get_script('vasp/subVASPatom_suncat.sh')

    @staticmethod
    def enmasse(elems     : L[int], curr : S[str],
                submitpth : str,
                beef      : str   = '',
                time      : int   = 10,
                sigma     : float = 0.01,
                sher      : bool  = True,
                ) -> None:
        """Submit many atomic jobs with the same parameters"""
        es = [Atomic(e,time,sigma) for e in elems]
        for e in es:
            e.submit(submitpth,sher=sher,curr=curr)

class Bulk(Job):
    """
    Multiple singlepoints at different strains to determine energy and bulkmod
    """
    def __init__(self,
                 bulk   : str,
                 time   : int   = 1,
                 sigma  : float = 0.01,
                 econv  : float = 1e-3,
                 strains: list  = list(range(-5,5)) # -9,10
                ) -> None:
        self.atoms   = read(bulk)
        self.name    = bulk[bulk.rfind('/')+1:bulk.find('.')]
        self.strains = strains
        calc = Calc(sigma = sigma, econv = econv, kpts = (10,10,10))
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Bulk<%s>'%self.name

    def submit(self,
               pth_  : str,
               sher  : bool ,
               retry : bool   = False,
               curr  : S[str] = set()
              ) -> None:
        pth = join(pth_,self.name)
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

                if sher: self._submit(spth)
                else:    self._submit_sunc(spth)

    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        incar = get_script('vasp/INCAR')
        incar_str = incar.format(Calc=self.calc)
        if self.calc.magmom is not None:
            incar_str+='\nMAGMOM = '+str(self.calc.magmom)
        with open(join(pth,'INCAR'),'w') as g:
            g.write(incar_str)
        kpts = get_script('vasp/KPOINTS')
        with open(join(pth,'KPOINTS'),'w') as g:
            g.write(kpts.format(*self.calc.kpts))

    @property
    def subsunc(self)->str:
        return get_script('vasp/subVASP_suncat.sh')

    @staticmethod
    def enmasse(pth    : str, submitpth : str, curr : S[str],
                beef   : str   = '',
                sher   : bool  = False,
                time   : int   = 1,
                sigma  : float = 0.01,
                econv  : float = 1e-3,
                strains: list  = list(range(-5,5)),
                retry  : bool  = False
                ) -> None:
        """ Create a bunch of bulks from a directory containing .traj files"""
        bulks = [Bulk(bulk=join(pth,p),time=time,sigma=sigma,econv=econv,strains=strains)
                    for p in listdir(pth) if p[-5:]=='.traj']
        for b in bulks:
            b.time = len(b.atoms)
            b.submit(pth_=submitpth, retry=retry, curr = curr, sher=sher)


#########################################
# PARSER #
def parse_elems(x : str) -> L[int]:
    if x == 'all': return allelems
    else:          return list(map(int, x.split()))

parser = ArgumentParser(description  = 'Submit some jobs',
                        allow_abbrev = True)

parser.add_argument('--time',
                    default = 1,
                    type    = int,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--sigma',
                    default = 0.01,
                    type    = float,
                    help    = 'Fermi temperature')

parser.add_argument('--econv',
                    default = 0.001,
                    type    = float,
                    help    = 'Energy convergence criterion ')

parser.add_argument('--src',
                    default = '',
                    type    = str,
                    help    = 'Path to bulk .traj files')

parser.add_argument('--target',
                    default = '',
                    type    = str,
                    help    = 'Path to where jobs will be submitted from')

parser.add_argument('--elems',
                    default = '',
                    type    = parse_elems,
                    help    = 'Either "all" or space separated list of positive integers')

parser.add_argument('--beef',
                    default = '',
                    help    = 'Copies a file into the working directory as BEEFoftheDay.txt')


def main()->None:
    '''Submit either bulk or element singlepoint calculations'''
    args = parser.parse_args()
    assert bool(args.src) ^ bool(args.elems), "Must be submiting Bulk or Atomic jobs"
    assert args.target,                       "Need a target location to submit jobs"

    calc = {attr:getattr(args,attr) for attr in ['time','sigma']}


    currstr = 'squeue -u ksb -o "%Z"' if sher else 'bjobs -o "job_name" -noheader'

    curr = set(check_output(currstr, encoding='UTF-8',shell=True).split())
    common = dict(submitpth=args.target,beef=args.beef,curr=curr)

    if args.elems:
        Atomic.enmasse(**{**common,**calc,**dict(sher=sher,elems=args.elems)})
    elif args.src:
        Bulk.enmasse(**{**common,**calc,**dict(sher=sher,pth=args.src)})
    else:
        raise ValueError("Must be submiting Bulk or Atomic jobs")

if __name__=='__main__':
    main()
