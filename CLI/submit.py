# External modules
from typing      import Union as U, List as L
from abc         import ABCMeta, abstractmethod
from os.path     import join,exists
from os          import chdir,system,mkdir,listdir,environ
from os          import mkdir
from argparse    import ArgumentParser

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
def get_script(s:str)->str:
    funcroot = environ['FUNCTIONALS_ROOT']
    with open(join(funcroot,'functionals/scripts/'+s),'r') as f:
        return f.read()

def safeMkdir(pth:str) -> None:
    if not exists(pth): mkdir(pth)

allelems = list(range(1,57)) + list(range(72,83))

class Calc(object):
    def __init__(self,
                 pw     : int   = 1000,
                 hund   : bool  = False,
                 david  : int   = 15,
                 sigma  : float = 0.01,
                 econv  : float = 5e-3,
                 dconv  : float = 1e-2,
                 nbands : int   = -8,
                 setups : str   = 'paw',
                 kpts   : U[int,float,tuple,list] = (1,1,1)
                ) -> None:

        if isinstance(kpts,(int,float)):
            self.kpts = {'density':kpts,'gamma':True} # type: U[dict,tuple]

        elif isinstance(kpts,(tuple,list)):
            self.kpts=tuple(kpts)
        else: raise ValueError

        self.pw = pw; self.hund = hund; self.david = david; self.sigma = sigma
        self.econv = econv; self.dconv = dconv; self.psp = setups;
        self.nbands = nbands

    def singlepoint(self,pth:str)->None:
        """ Write a singlepoint calculation to a file """
        singlepoint = get_script('singlepoint.py')
        with open(join(pth,'gpawrun.py'),'w') as g:
            g.write(singlepoint.format(Calc=self))

class Job(metaclass = ABCMeta):
    def __init__(self, calc : Calc, time : int = 1) -> None:
        assert time > 0
        self.calc = calc
        self.time = time

    def __repr__(self)->str:
        return str(self)

    def _submit(self, pth : str) -> None:
        self.calc.singlepoint(pth)
        bash = join(pth,'subGPAW.sh')
        sub  = get_script('subGPAW.sh')
        with open(bash,'w') as g:
            g.write(sub.format(Job=self))
        chdir(pth);
        system('chmod 755 '+bash)
        system('sbatch '+bash)

    @abstractmethod
    def submit(self,pth:str)->None:
        raise NotImplementedError

class Atomic(Job):
    """
    Calculation of the energy for an isolated atom
    """
    def __init__(self,
                 elem   : int,
                 time   : int   = 10,
                 sigma  : float = 0.01,
                 econv  : float = 1e-3,
                 dconv  : float = 1e-3,
                ) -> None:

        self.elem = elem
        calc = Calc(sigma=sigma, econv=econv, dconv=dconv, hund=True, kpts=(1,1,1))
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Atoms<%d>'%(chemical_symbols[self.elem])

    def submit(self, pth_ : str, retry : bool = False) -> None:
        pth = join(pth_,chemical_symbols[self.elem])
        safeMkdir(pth)
        if retry or not exists(join(pth,'xccontribs.txt')):

            atoms = Atoms(numbers    = [self.elem],
                          cell        = [[10.1,0.2,0.1],
                                         [0.2,10.2,0.01],
                                         [0.5,0.1,10.3]],
                          positions    = [[0.01,0.02,0.03]])
            write(join(pth,'init.traj'),atoms)
            self._submit(pth)

    @staticmethod
    def enmasse(elems     : L[int],
                submitpth : str,
                beef      : str   = '',
                time      : int   = 10,
                sigma     : float = 0.01,
                econv     : float = 1e-3,
                dconv     : float = 1e-3,
                ) -> None:
        """Submit many atomic jobs with the same parameters"""
        es = [Atomic(e,time,sigma,econv,dconv) for e in elems]
        for e in es:
            e.submit(submitpth)

class Bulk(Job):
    """
    Multiple singlepoints at different strains to determine energy and bulkmod
    """
    def __init__(self,
                 bulk   : str,
                 time   : int   = 1,
                 sigma  : float = 0.01,
                 econv  : float = 1e-3,
                 dconv  : float = 1e-3,
                 strains: list  = list(range(-9,10))
                ) -> None:
        self.atoms   = read(bulk)
        self.name    = bulk[bulk.rfind('/')+1:bulk.find('.')]
        self.strains = strains
        calc = Calc(sigma = sigma, econv = econv, dconv = dconv, kpts = 5)
        super().__init__(calc=calc,time=time)

    def __str__(self)->str:
        return 'Bulk<%s>'%self.name

    def submit(self, pth_ : str, retry : bool = False) -> None:
        pth = join(pth_,self.name)
        safeMkdir(pth)
        opt_cell = self.atoms.get_cell()
        for strain in self.strains:
            cell = opt_cell * (1 + strain/100.)
            spth = join(pth,'strain_%d'%strain)
            if retry or not exists(join(spth,'xccontribs.txt')):
                safeMkdir(spth)
                strained_atoms = self.atoms.copy()
                strained_atoms.set_cell(cell, scale_atoms = True)
                write(join(spth,'init.traj'),strained_atoms)
                self._submit(spth)

    @staticmethod
    def enmasse(pth    : str, submitpth : str,
                beef   : str   = '',
                time   : int   = 1,
                sigma  : float = 0.01,
                econv  : float = 1e-3,
                dconv  : float = 1e-3,
                strains: list  = list(range(-10,11)),
                retry  : bool  = False
                ) -> None:
        """ Create a bunch of bulks from a directory containing .traj files"""
        bulks = [Bulk(bulk=join(pth,p),time=time,sigma=sigma,econv=econv,dconv=dconv,strains=strains)
                    for p in listdir(pth) if p[-5:]=='.traj']
        for b in bulks:
            b.submit(submitpth, retry)


#########################################
# PARSER #
def parse_elems(x : str) -> L[int]:
    if x == 'all': return allelems
    else:          return list(map(int, x.split()))

parser = ArgumentParser(description  = 'Submit some jobs',
                        allow_abbrev = True)

parser.add_argument('--time',
                    default = 5,
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

parser.add_argument('--dconv',
                    default = 0.001,
                    type    = float,
                    help    = 'Density convergence criterion')

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


    calc = {attr:getattr(args,attr) for attr in ['time','sigma','econv','dconv']}

    common = {'submitpth':args.target,'beef':args.beef}

    if args.elems:
        Atomic.enmasse(elems=args.elems,**common,**calc) # type: ignore
    elif args.src:
        Bulk.enmasse(pth=args.src,**common,**calc) # type: ignore
    else:
        raise ValueError("Must be submiting Bulk or Atomic jobs")

if __name__=='__main__':
    main()
