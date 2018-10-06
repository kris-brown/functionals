from sys import version_info
if version_info[0] == 3:
    from typing import Any
import sys
from json 	     import dump
from gpaw 		 import GPAW, PW, FermiDirac, restart, Davidson, Mixer, MixerSum, MixerDif # type: ignore
from gpaw.xc.bee import BEEFEnsemble # type: ignore
from ase 	     import Atoms        # type: ignore
from ase.io      import read         # type: ignore

"""
To format file, read this into string variable and write var.format(Calc)

Files needed:
    init.traj
Files produced:
    inp.gpw
    energy.txt
    xccontribs.txt
"""

Calc = None # type: Any

atoms = read('init.traj')

atoms.calc = GPAW(setups       = '{Calc.psp}'
                 ,mode         = PW({Calc.pw})
                 ,xc           = 'PBE'
                 ,hund         = {Calc.hund}
                 ,eigensolver  = Davidson({Calc.david})
                 ,kpts         = {Calc.kpts}
                 ,occupations  = FermiDirac({Calc.sigma},fixmagmom={Calc.hund})
                 ,nbands       = {Calc.nbands}
                 ,txt          = 'pbe.txt')

atoms.get_potential_energy()

atoms.calc.write('inp.gpw', mode='all')

atoms,calc = restart('inp.gpw')

calc.set(xc           = 'mBEEF'
		,setups       = '{Calc.psp}'
        ,mode         = PW({Calc.pw})
        ,hund         = {Calc.hund}
        ,kpts         = {Calc.kpts}
        ,eigensolver  = Davidson({Calc.david})
        ,occupations  = FermiDirac({Calc.sigma},fixmagmom={Calc.hund})
        ,nbands       = {Calc.nbands}
        ,convergence  = {{'energy': {Calc.econv},'density':{Calc.dconv},'eigenstates':1e80}}
        ,txt          = 'mbeef.txt')

e = atoms.get_potential_energy()

with open('energy.txt','w') as f:
    f.write(str(e))

beef = BEEFEnsemble(calc)

xc = beef.mbeef_exchange_energy_contribs()

with open('xccontribs.txt','w') as f:
    dump(xc.tolist(),f)
