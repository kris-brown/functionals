import os
from ase.io import read,write  # type: ignore
from ase.eos  import EquationOfState # type: ignore
import matplotlib # type: ignore
matplotlib.use('Qt4Agg')
from pymatgen.io.ase 			import AseAtomsAdaptor 	  # type: ignore
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer # type: ignore

"""
Optimize a cell lattice by straining it (i.e. assuming only one free DoF)
"""
############################################################################################



#######
# INPUTS
######
pw       = 1000
exclude  = []  # type: list
lower    = -7
upper    = 7

######

jobstr = """
from gpaw import GPAW,PW,Davidson,FermiDirac,restart
from ase.io import read
from gpaw.xc.bee import BEEFEnsemble # type: ignore
from json import dump
#####################
atoms = read('init.traj')
mdict = {{'Ni':0.8,'Co':1.9,'Fe':2.4}}
atoms.set_initial_magnetic_moments([mdict.get(x,0) for x in atoms.get_chemical_symbols()])
atoms.calc = GPAW(setups       = 'paw'
                 ,mode         = PW( {0} )
                 ,xc           = 'PBE'
                 ,eigensolver  = Davidson(15)
                 ,kpts         = {{'density': 4.5, 'gamma': True}}
                 ,occupations  = FermiDirac(0.01)
                 ,nbands       = -8
                 ,convergence  = {{'energy': 0.05,'density':1e-2,'eigenstates':1e80}}
                 ,txt          = 'pbe.txt')
atoms.get_potential_energy()
atoms.calc.write('inp.gpw', mode='all')
atoms,calc = restart('inp.gpw')
calc.set(xc='mBEEF'
        ,mode         = PW({0})
        ,eigensolver  = Davidson(15)
        ,kpts         ={{'density': 4.5, 'gamma': True}}
        ,nbands       = -8
        ,occupations     = FermiDirac(0.001)
        ,convergence     = {{'energy': 0.005,'density':1e-3,'eigenstates':1e80}}
        ,txt             = 'mbeef.txt')
e = atoms.get_potential_energy()
with open('energy.txt','w') as f:
    f.write(str(e))
beef = BEEFEnsemble(calc)
xc = beef.mbeef_exchange_energy_contribs()
with open('xccontribs.txt','w') as f:
    dump(xc.tolist(),f)

""".format(pw)

metastr = """#!/bin/bash
#SBATCH -p iric,suncat
#SBATCH --job-name=myjob
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --time=00:55:00
#SBATCH --nodes=1
#SBATCH --mincpus=16
#SBATCH --mem-per-cpu=4000
#SBATCH  --mail-user=ksb@stanford.edu

NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`

source /scratch/users/ksb/gpaw/paths.bash
export GPAW_SETUP_PATH=/scratch/users/ksb/atomization/gpaw-setups-0.8.7929

mpirun -n 16 gpaw-python opt.py
"""
#############
# Main Script
#----------
def main()->None:
    atoms    = read('init.traj')
    pmg      = AseAtomsAdaptor().get_structure(atoms)
    sg       = SpacegroupAnalyzer(pmg)
    pmgunit  = sg.get_conventional_standard_structure()
    atoms    = AseAtomsAdaptor().get_atoms(pmgunit)
    opt_cell = atoms.get_cell()
    rootdir  = os.getcwd()
    vol      = []
    done     = True
    strains  = [x for x in range(lower,upper) if x not in exclude]
    for strain in strains:
        strained_atoms = atoms.copy()
        strained_atoms.set_cell(opt_cell*(1+strain/100.),scale_atoms=True)
        vol.append(strained_atoms.get_volume())
        newdir = os.getcwd()+'/strain_%d'%strain

        try:
            os.mkdir(newdir)
        except FileExistsError:
            pass

        os.chdir(newdir)

        write('init.traj',strained_atoms)
        with open('opt.py','w') as f:
            f.write(jobstr)
        with open('sub.sh','w') as f:
            f.write(metastr)

        def is_real(x:str)->bool:
            return os.path.exists(x) and os.path.getsize(x)>100

        if not is_real('xccontribs.txt'):
            print('sbatch in '+newdir)
            os.system('sbatch -J $PWD sub.sh')
            done = False
        os.chdir(rootdir)


    if not done:
        print('\tStill had jobs to submit, analyze lattice opt results another day')
        return None
    ################################
    #Analyzing dE/dV curve
    #-------------------------------
    eng = []
    for strain in strains:
        with open('strain_%d/energy.txt'%strain,'r') as f:
            eng.append(float(f.read()))
    eos = EquationOfState(vol,eng)
    v0, e0, b = eos.fit()
    print('Volume %f, Energy %f, Bulkmod %f'%(v0,e0,b))
    eos.plot()

if __name__=='__main__':
    main()
