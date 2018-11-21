# this script is meant to be called with sherlock's gpaw-python
# which is a python2 implementation
from sys     import argv
from os.path import exists
from json    import load, dump
from numpy   import array,gradient,tile,pi,vectorize,sum,average  # type: ignore
from numpy.linalg import norm # type: ignore
from numpy.random import choice # type: ignore

'''
About the awkward structure of this file
#-----------------------------------------
Import statements for modules that can only be run on either Sherlock or
local machine (but not both) must be hidden below in the code.

    - It is currently impossible to generate the s + alpha distribution on a local
        machine (because GPAW cannot be installed).

    - It is currently impossible to visualize the data on sherlock (cannot pip install
        anything, including plotly/matplotlib, due to unknown problems. Furthermore, unknown
        problems with X11 forwarding prevent visualization even when matplotlib is there)

Thus half of the work is done on Sherlock (generating a JSON file, which gets
  scp'd to local machine), and half is done on local machine:
    - When this script is called with one argument (or none), the first half
        gets executed.
    - With two args, the second half is executed.


What kind of results should be shown?
#--------------------------------------

The average s or alpha over volume is pretty meaningless:
    - at very low densities, s blows up due inversely scaling with r^(4/3)
    - s essentially becomes proportional to the amount of vacuum: not good

Weighting by energy seems more useful, but the resulting product of these
    numbers is hard to interpret.

Instead, I randomly sample s or a weighted by density to make a histogram
'''
##############################################################################

if __name__=='__main__':

    r3 = range(3) # shorthand

    # Visualize Results
    #------------------
    if len(argv) > 2:
        from plotly            import tools # type: ignore
        from plotly.graph_objs import Histogram as Hist,Layout,Figure # type: ignore
        from plotly.offline    import plot # type: ignore

        # Load input
        #-----------
        assert exists(argv[1]), 'Point to a JSON output from the other part of this code'

        with open(argv[1],'r') as f:
            s,a,pden = map(array,load(f))

        # Randomly sample with weighting
        #-------------------------------
        ss  = [choice(a=s.flatten(),p=pden) for _ in range(1000)]
        as_ = [choice(a=a.flatten(),p=pden) for _ in range(1000)]

        # Plotly minutiae
        #----------------
        fig = tools.make_subplots(rows=1, cols=2, specs=[[{}, {}]],
                                  shared_xaxes=False, shared_yaxes=False,
                                  )

        args   = dict(histnorm ='probability',nbinsx = 500) # common kwargs
        data   = [Hist(x=ss,name='s',**args),Hist(x=as_,name='alpha',xaxis='x2',**args)]
        fig.append_trace(data[0], 1, 1)
        fig.append_trace(data[1], 1, 2)

        layout = Layout(xaxis=dict(title='s'),xaxis2=dict(title='a'))
        plot(fig)

    else:

        from gpaw    import restart  # type: ignore

        # Get a gpw file
        #---------------

        if len(argv) > 1:
            assert exists(argv[1]), 'Need path to .gpw file!'
            pth = argv[1]
        else:
            pth = '/scratch/users/ksb/functionals/data/example.gpw' # Cubic aluminum, strained geometry

        # Extract atoms and density from file
        #-----------------------------------
        atoms,calc = restart(pth)
        den        = calc.get_pseudo_valence_density() # 3D numpy tensor
        den[den<0] = 1e-60 # sometimes density is (very slightly) negative, wow.
        cell       = atoms.get_cell()

        # Determine grid spacing for each coordinate
        #-------------------------------------------
        n,ny,nz  = den.shape
        vecs     = [norm(cell[i]) for i in r3] # cell vector lengths
        spacing  = [vec/gpoints for vec,gpoints in zip(vecs,den.shape)]

        # Make replicated cell in all three dimensions to create effective P.B.C
        # for the lone cube (out of 27) that is in the center
        #-----------------------------------------------------------------------
        d3 = tile(den,(3,3,3))

        # Compute gradients with constant spacing
        # (though potentially different along each axis)
        #-----------------------------------------------
        d1 = gradient(d3,*spacing) # list of three 3Nx3Nx3N tensors
        d2 = [gradient(d,spacing[ax],axis=ax) for d in d1 for ax in r3] # nine 3Nx3Nx3N tensors


        # Compute Euclidean norm of all tensors, then slice out the center cube
        #---------------------------------------------------------------------------
        slicer = [slice(n, 2*n, 1)] * 3 # center cube in a 3Nx3Nx3N tensor
        grad   = norm([d[slicer] for d in d1],axis=0) # |1st derivative|
        grad2  = norm([d[slicer] for d in d2],axis=0) # |2nd derivative|

        # Pointwise computation of s and alpha (r == rho)
        #-----------------------------------------------
        # Helpers
        pi2   = pi ** 2.
        kf    = lambda r:    (3. * pi2 * r)**(1./3.)
        t_ueg = lambda r:     0.3 * (3. * pi2)**(2./3.) * r**(5./3.)
        t_w   = lambda r, dr: dr**2 / 8 / r
        t     = lambda d2r:   d2r  # ???
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
