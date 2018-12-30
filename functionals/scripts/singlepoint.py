from sys import version_info
if version_info[0] == 3:
    from typing import Any

from json        import dump
from gpaw        import GPAW, PW, FermiDirac, restart, Davidson, Mixer # type: ignore ##, MixerSum, MixerDif
from gpaw.xc.bee import BEEFEnsemble # type: ignore
from ase.io      import read         # type: ignore
from numpy        import array,gradient,tile,pi,vectorize,sum,average  # type: ignore
from numpy.linalg import norm # type: ignore
from numpy.random import choice # type: ignore

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


    atoms.calc.write('inp.gpw', mode='all')

    atoms,calc = restart('inp.gpw')

#### Record Energy ####
e = atoms.get_potential_energy()
with open('energy.txt','w') as f:
    f.write(str(e))

#### Record basis function contribs ####

beef = BEEFEnsemble(calc)

xcc = beef.mbeef_exchange_energy_contribs()

with open('xccontribs.txt','w') as f:
    dump(xcc.tolist(),f)


#### Compute s and alpha distributions ####

r3         = range(3) # shorthand
wf         = calc.get_pseudo_wave_function()   # 3D numpy tensor
den        = calc.get_pseudo_valence_density() # 3D numpy tensor
den[den<0] = 1e-60 # sometimes density is (very slightly) negative, wow.
cell       = atoms.get_cell()
vecs       = [norm(cell[i]) for i in r3] # cell vector lengths


# Get gradients
#--------------
gradients  = [None,None] # fill this in the for loop

for counter,tensor in enumerate([wf,den]):
    # Determine grid spacing for each coordinate
    #-------------------------------------------
    spacing  = [vec/gpoints for vec,gpoints in zip(vecs,tensor.shape)]

    # Make replicated cell in all three dimensions to create effective P.B.C
    # for the lone cube (out of 27) that is in the center
    #-----------------------------------------------------------------------
    bigcube = tile(tensor,(3,3,3))                      # 3N x 3N x 3N tensor

    # trim so that there is just a border of thickness 2
    #-----------------------------------------------------------------------
    slicer1 = [slice(n-2, 2*n+2) for n in tensor.shape] # 3N -> N+4 - cube tensor
    slicer2 = [slice(2, -2)] * 3                        # N+4 -> N - cube tensor

    trimmed = bigcube[slicer1]  # do gradient on as small of tensor as poss.

    # Compute gradients with constant spacing
    # (though potentially different along each axis)
    #-----------------------------------------------
    d1 = gradient(trimmed,*spacing) # LIST of three N+4 cube tensors

    # Trim off edges and compute Euclidean norm of all tensors,
    #---------------------------------------------------------------------------
    gradients[counter]  = norm([d[slicer2] for d in d1],axis=0) # |Derivative|

grad2,grad = gradients

# Pointwise computation of s and alpha (r == rho)
#-----------------------------------------------

# Helpers
pi2   = pi ** 2.
kf    = lambda r:     (3. * pi2 * r)**(1./3.)
t_ueg = lambda r:     0.3 * (3. * pi2)**(2./3.) * r**(5./3.)
t_w   = lambda r, dr: dr**2 / 8 / r
t     = lambda d2r:   0.5* d2r**2

# Important functions
s     = lambda r, dr:       dr/(2. * kf(r) * r)              # type: ignore
alpha = lambda r, dr, d2r: (t(d2r) - t_w(r,dr))/t_ueg(r)     # type: ignore

# Compute the pointwise operations operate over matrices
#-------------------------------------------------------
s_mat = vectorize(s)(den,grad)           # NxNxN tensor
a_mat = vectorize(alpha)(den,grad,grad2) # NxNxN tensor
pden  = den / sum(den)

rpt = lambda x: [round(x(),2) for _ in range(1000)]
dists = [rpt(lambda:choice(a=mat.flatten(),p=pden.flatten())) for mat in [s_mat,a_mat]] # type: ignore

with open('distribs.json','w') as f:
    dump(dists,f)
