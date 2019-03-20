from typing         import Tuple as T, Dict as D, List as L, Any
from ast            import literal_eval
from numpy          import array,vstack,zeros,ones,concatenate,empty,sum,average,abs,split # type: ignore
from numpy.linalg   import inv   # type: ignore
'''
Data Preprocessing for nonlinear fitting
'''
T6 = T[array,array,array,array,array,array]
################################################################################

def weight(w_bm:float, w_lat:float, arrays:tuple) -> T[array,array]:
    x_ce,x_bm,x_lat,y_ce,y_bm,y_lat = arrays
    w_bm_ = w_bm * average(abs(y_bm))/average(abs(y_ce))
    w_lat_ = w_lat * average(abs(y_lat))/average(abs(y_ce))
    return (vstack((x_ce, w_bm_*x_bm, w_lat_*x_lat)),
            concatenate((y_ce, w_bm_*y_bm, w_lat_*y_lat)))

def process_data(data : L[D[str,Any]]) -> T6:
    '''
    Take information about a particular chemical system and return useful
    data structures for fitting based on cohesive energy, lattice constant, and
    bulk modulus
    '''
    if data:

        processed   = [_extract(d) for d in data]

        import pdb;pdb.set_trace()

        cohesiveX,cohesiveY,bmX,bmY,latticeX,latticeY = zip(*processed)
        X0,X1,X2  = map(vstack,[cohesiveX,bmX,latticeX])
        Y0,Y1,Y2  = map(array,[cohesiveY,bmY,latticeY])
    else:
        X0,X1,X2 = zeros((1,8**2)),empty((0,8**2)),empty((0,8**2))
        Y0,Y1,Y2 = [0],empty((0)),empty((0))

    return X0,X1,X2,Y0,Y1,Y2

def _extract(d : dict) -> T6:
    '''
    Calls all of the preprocessing scripts
    '''
    cohesiveX,cohesiveY = _get_cohesive_x_target(d)
    bmX, bmY, latticeX, latticeY = _get_bm_lattice_x_target(d)

    return cohesiveX,cohesiveY,bmX,bmY,latticeX,latticeY

def _get_cohesive_x_target(d : dict)->T[array,array]:
    '''
    Convert a row data into a row of coefficients and a target for cohesive data
    fitting. Requires data to have the following keys:
    - composition, atomic_contribs, atomic_energies, coefs, bulk_ratio,
      bulk_energy, bulk_contribs, expt_cohesive_energy

    returns a 5x64 matrix and a length-5 vector
    '''

    # Extract info from dictionary
    comp   = literal_eval(d['composition'])             # type: D[int,int]
    raw_ac = literal_eval(d['atomic_contribs']).items() # int -> 64 element array
    raw_ae = literal_eval(d['atomic_energies']).items() # int -> float

    coefs  = array(literal_eval(d['coefs']))            # 64 element array
    ratio  = int(d['bulk_ratio'])                       # actual system : reference system
    ex_ce  = float(d['expt_cohesive_energy'])           # experimental E atom - E bulk
    e_bulk = float(d['bulk_energy']) / ratio            # calculated energy of bulk reference system

    contribs = array(literal_eval(d['bulk_contribs'])).flatten()

    x_bulks = array(split(contribs / ratio,5))                           # 5 1x64 arrays

    # Analysis
    #---------
    # Get coefficients representing the change in exchange contributions
    atom_contribs = {i:array(xs) for i,xs in raw_ac} # int -> 8x8 matrix
    atom_energies = {i:float(e)  for i,e  in raw_ae} # int -> float

    e_atom  = sum([atom_energies[e]*num for e,num in comp.items()]) # float
    x_atom  = sum([atom_contribs[e]*num for e,num in comp.items()],axis=0) # 8x8 matrix
    x_atoms = split(x_atom.flatten(),5) # 5 1x64 arrays

    dx = (x_atoms - x_bulks) # 5 x 64 element array, THIS IS WHAT WE ARE FITTING



    # Get target to fit the above to: JUST the ex_component of cohesive energy
    ex_atoms = [coefs @ x for x in x_atoms] # Use the BEEF coefficients
    ex_bulks = [coefs @ x for x in x_bulks] # from this particular calculator

    nonx_e_atom = array([e_atom - ex for ex in ex_atoms])
    nonx_e_bulk = array([e_bulk - ex for ex in ex_bulks])
    target      = ex_ce  - (nonx_e_atom - nonx_e_bulk) # Just Ex_atom - Ex_bulk

    return dx,target

def _get_bm_lattice_x_target(d : dict) -> T[array,array,array,array]:
    '''
    Fit energies to quadratic form using linear algebra
    ---------------------------------------------------
    Solve b = A x
        with: x = (Aᵗ·A)⁻¹ · Aᵗ · b
            - x = vector of three elements: x2,x1,x0
                - such that: Energy = x2*Vol²+x1*Vol+x0
            - b = contribs·coefs + e_nonx
                - this is a function of the vector being fit
            - A = a "vandermonde" matrix generated purely by the volume data

    so x = (Aᵗ·A)⁻¹·Aᵗ·(contribs·coefs + e_nonx)
        or x[i](coefs) = VECTOR · COEFS + CONST
            - VECTOR = i'th row of (Aᵗ·A)⁻¹·Aᵗ·contribs
            - CONST  = i'th element of 1 x 3 vector: (Aᵗ·A)⁻¹·Aᵗ·contribs · e_nonx

    Curvature prediction error = (x[2](coefs) - expt_curv)
    Lattice prediction error   = (x[1](coefs) / 2*expt_curv) - expt_volume
        - note that minimum of parabola is at -x1/2*x2

    In order to make the error some A·x - b, we need:
        Curvature:
            - A = x[2] VECTOR
            - b = expt_curv - x[2] CONST
        Lattice:
            - A = x[1] VECTOR / 2*expt_curv
            - b = expt_vol  - (x[1] CONST / 2*expt_curv)

    We need this for all 5 functionals, so our return types are
    A_bm --- 5 x 64
    '''

    # Common stuff to preprocessing both BM and Lattice data
    #----------------------------------------------------------

    import pdb;pdb.set_trace()

    energies = array(d['energy_vector'])   # 5 element array
    volumes  = array(d['volume_vector'])    # 5 element array

    contribs = array(d['contrib_vector']).reshape((5,5,64)) # 5 x 5 x 64

    bm          = d['expt_bm']*10**9     # experimental value, Pa or N / m²
    expt_volume = d['expt_volume']*(10**-30)       # experimental volume, m^3

    import pdb;pdb.set_trace()

    curv_vecs, curv_consts, lat_vecs, lat_consts = [],[],[],[]

    for c in contribs:
        # Get experimental curvature to E vs V
        #--------------------------------------
        curv_      = bm / expt_volume                   # experimental d²E/dV², J/m^6 = N/m^5
        expt_curv  = curv_ * (10**-60) * (6.242*10**18) # eV / A^6

        # allows us to convert xs -> energy
        coefs = array(literal_eval(d['coefs']))

        e_nonx  = energies - c @ coefs # len-5 vector

        # also because we only need "a", we only dot the last row with the coef (col) vec

        vander = vstack((ones(len(volumes)),volumes,volumes**2)).T  # 5 x 3
        vinv   = inv(vander.T @ vander)                             # 3 x 3
        solver =  vinv @ vander.T                                   # 3 x 5
        vecs    = solver @ c                                        # 3 x 64*5
        constvec  = solver @ e_nonx                                 # 1 x 3

        curv_vec   = vecs[2],
        curv_const = expt_curv - constvec[2],
        lat_vec    = vecs[1] / 2*expt_curv
        lat_const  = expt_volume -  constvec[1] / 2*expt_curv

        import pdb;pdb.set_trace()
        curv_vecs.append(curv_vec)
        curv_consts.append(curv_const[0])
        lat_vecs.append(lat_vec)
        lat_consts.append(lat_const)

    return  curv_vec, curv_consts, lat_vecs, lat_consts
