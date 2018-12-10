# this script is meant to be called with sherlock's gpaw-python
# which is a python2 implementation
from sys     import argv
from os.path import exists
from json    import dump
from numpy   import array,gradient,tile,pi,vectorize,sum,average  # type: ignore
from numpy.linalg import norm # type: ignore

'''
'''
##############################################################################

if __name__=='__main__':
    from gpaw import restart  # type: ignore

    r3 = range(3) # shorthand

    # Get a gpw file
    #---------------

    if len(argv) > 1:
        assert exists(argv[1]), 'Need path to .gpw file!'
        pth = argv[1]
    else:
        pth = '/scratch/users/ksb/functionals/data/example.gpw' # Cubic aluminum, strained geometry

    # Extract atoms, wf, and density from file
    #---------------------------------------
    atoms,calc = restart(pth)

    print('fix treatment of pseudo_wave_function')
    '''
    T = -1/2  < psi | del^2 | psi>
      integrate by parts
    T = | del psi | ^2
    '''
    # look at lowest energy band
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

    with open('out.json','w') as f:
        dump([x.flatten().tolist() for x in [s_mat,a_mat,pden]],f)
