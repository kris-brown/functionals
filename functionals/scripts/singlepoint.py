from sys import version_info
if version_info[0] == 3:
    from typing import Any
import sys
from json          import dump
from gpaw          import GPAW, PW, FermiDirac, restart, Davidson, Mixer, MixerSum, MixerDif # type: ignore
from gpaw.xc.bee import BEEFEnsemble # type: ignore
from ase          import Atoms        # type: ignore
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

common_calc_params = {{
    'setups'      : '{Calc.psp}',
    'hund'        : {Calc.hund},
    'mode'        : PW({Calc.pw}),
    'kpts'        : {Calc.kpts},
    'maxiter'     : 1000,
    'mixer'       : Mixer(0.01,  5,  100),
    'eigensolver' : Davidson({Calc.david}),
    'nbands'      : {Calc.nbands},
    'occupations' : FermiDirac({Calc.sigma},fixmagmom = {Calc.hund}),
    'convergence' : {{'energy': {Calc.econv}, 'density':{Calc.dconv}, 'eigenstates':1e80}},
}}

xcs = ['PBE','mBEEF']

for xc in xcs:
    atoms.calc = GPAW(xc  = xc,
                      txt = xc+'.txt',
                      **common_calc_params)

    e = atoms.get_potential_energy()

    atoms.calc.write('inp.gpw', mode='all')

    atoms,calc = restart('inp.gpw')

with open('energy.txt','w') as f:
    f.write(str(e))

beef = BEEFEnsemble(calc)

xcc = beef.mbeef_exchange_energy_contribs()

with open('xccontribs.txt','w') as f:
    dump(xcc.tolist(),f)
