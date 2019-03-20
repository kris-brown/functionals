from typing     import Any, List as L, Dict as D, Tuple as T
from argparse   import ArgumentParser
from os.path    import join, exists
from itertools  import product as prod
from os         import environ, system, listdir
from shutil     import copyfile
from pathlib    import Path
from random     import choice
from ast        import literal_eval
from json       import load,dump, loads
from numpy      import array,inf,sum,vstack,ones # type: ignore
from numpy.linalg import inv   # type: ignore
from MySQLdb    import connect,Connection # type: ignore

from functionals.fit.fit import Fit,sqlselect

'''
Submit fitting jobs, provided a DB connection to the DB with DFT jobs
'''
###############################################################################
# Constants
###########

# Queries
#------------
q1 = 'SELECT fitparams_id  FROM fitparams'
q2 = 'SELECT DISTINCT calc FROM expt'

def main(db_ : str, pth : str) -> None:
    '''
    For every completed "expt" and every single set of fitting parameters,
    set up 5 fitting jobs
    '''

    with open(db_,'r') as f: conn = connect(**load(f),  autocommit  = True)

    params,calcs = [sqlselect(conn,x) for x in [q1,q2]]

    for (fp,),(calc,) in prod(params,calcs):
        for decay in range(5):
            fit = Fit.from_db(db=db_,fp_id=fp,calc_id=calc,decay=decay)
            root = join(pth,fit.uid()[:10],str(decay))
            Path(root).mkdir(parents=True, exist_ok=True)
            fit.write(root)

def sub(pth : str, time : int, retry : bool, local : bool) -> None:
    dirs = listdir(pth)
    for d in dirs:
        for dd in listdir(join(pth,d)):
            dir   = join(pth,d,dd)
            if retry or not exists(join(dir,'result.json')):
                sunc  = choice(['','','','2','2','3'])
                act   = 'python runfit.py' if local else 'bsub -n 1 -W{}:09 -q suncat{} subfit.sh'.format(time,sunc)
                cmd   = 'cd {}; '.format(dir) + act
                system(cmd)


# Parser
########
parser = ArgumentParser(description  = 'Submit some fitting jobs',
                        allow_abbrev = True)

parser.add_argument('--db', type    = str,
                    default = '/Users/ksb/Documents/JSON/functionals.json',
                    help    = 'Path to JSON with DB connection info')

parser.add_argument('--pth', type    = str,
                    default = '/Users/ksb/scp_tmp/fit',
                    help    = 'Path to where fitting jobs will be performed')

parser.add_argument('--sub',
                    help='Add anything to do submit command instead')

parser.add_argument('--time', type    = int,
                    default = 1,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--retry', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--local', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.sub:
        sub(args.pth,args.time,args.retry,args.local)
    else:
        main(args.db,args.pth)



# def cohesive(comp_   : str,
#              raw_ac_ : D[int,array],
#              raw_ae_ : str,
#              coefs   : array,
#              ratio   : int,
#              tar     : float,
#              ebulk_  : float,
#              bcontribs_: array,
#             ) -> T[list,list]:
#     '''
#     Convert a row data into a row of coefficients and a target for cohesive data
#     fitting. Requires data to have the following keys:
#     - composition, atomic_contribs, atomic_energies, coefs, bulk_ratio,
#       bulk_energy, bulk_contribs, expt_cohesive_energy
#
#     returns a 5x64 matrix and a length-5 vector
#     '''
#
#     # Extract info from dictionary
#     comp   = literal_eval(comp_)             # type: D[int,int]
#     raw_ac = raw_ac_.items() # int -> 64 element array
#     raw_ae = literal_eval(raw_ae_).items() # int -> float
#
#     ex_ce  = float(tar)           # experimental E atom - E bulk
#     e_bulk = float(ebulk_) / ratio            # calculated energy of bulk reference system
#
#     contribs = bcontribs_.flatten()
#
#     x_bulk = contribs / ratio # 64 element array
#
#     # Analysis
#     #---------
#     # Get coefficients representing the change in exchange contributions
#     atom_contribs = {i:array(xs) for i,xs in raw_ac} # int -> 8x8 matrix
#     atom_energies = {i:float(e)  for i,e  in raw_ae} # int -> float
#
#     e_atom  = sum([atom_energies[e]*num for e,num in comp.items()]) # float
#     x_atom  = sum([atom_contribs[e]*num for e,num in comp.items()],axis=0) # 8x8 matrix
#
#     dx = (x_atom - x_bulk) # 64 element array, THIS IS WHAT WE ARE FITTING
#
#     # Get target to fit the above to: JUST the ex_component of cohesive energy
#     ex_atom = coefs @ x_atom # Use the BEEF coefficients
#     ex_bulk = coefs @ x_bulk # from this particular calculator
#
#     nonx_e_atom = e_atom - ex_atom
#     nonx_e_bulk = e_bulk - ex_bulk
#     target      = ex_ce  - (nonx_e_atom - nonx_e_bulk) # Just Ex_atom - Ex_bulk
#
#
#     return dx.tolist(),target
#
# def bm_lat(engs     : str,
#            vols     : str,
#            contribs : list,
#            expt_bm  : float,
#            expt_vol : float,
#            coefs    : array,
#           ) -> T[array,array,array,array]:
#     '''
#     Fit energies to quadratic form using linear algebra
#     ---------------------------------------------------
#     Solve b = A x
#         with: x = (Aᵗ·A)⁻¹ · Aᵗ · b
#             - x = vector of three elements: x2,x1,x0
#                 - such that: Energy = x2*Vol²+x1*Vol+x0
#             - b = contribs·coefs + e_nonx
#                 - this is a function of the vector being fit
#             - A = a "vandermonde" matrix generated purely by the volume data
#
#     so x = (Aᵗ·A)⁻¹·Aᵗ·(contribs·coefs + e_nonx)
#         or x[i](coefs) = VECTOR · COEFS + CONST
#             - VECTOR = i'th row of (Aᵗ·A)⁻¹·Aᵗ·contribs
#             - CONST  = i'th element of 1 x 3 vector: (Aᵗ·A)⁻¹·Aᵗ·contribs · e_nonx
#
#     Curvature prediction error = (x[2](coefs) - expt_curv)
#     Lattice prediction error   = (x[1](coefs) / 2*expt_curv) - expt_volume
#         - note that minimum of parabola is at -x1/2*x2
#
#     In order to make the error some A·x - b, we need:
#         Curvature:
#             - A = x[2] VECTOR
#             - b = expt_curv - x[2] CONST
#         Lattice:
#             - A = x[1] VECTOR / 2*expt_curv
#             - b = expt_vol  - (x[1] CONST / 2*expt_curv)
#
#     We need this for all 5 functionals, so our return types are
#     A_bm --- 5 x 64
#     '''
#
#     # Common stuff to preprocessing both BM and Lattice data
#     #----------------------------------------------------------
#
#     energies = array(loads(engs))   # 5 element array
#     volumes  = array(loads(vols))   # 5 element array
#
#     contribs = array(contribs).reshape((5,64)) # 5 x 64
#
#     bm          = expt_bm*10**9     # experimental value, Pa or N / m²
#     expt_volume = expt_vol*(10**-30)       # experimental volume, m^3
#
#     # Get experimental curvature to E vs V
#     #--------------------------------------
#     curv_      = bm / expt_volume                   # experimental d²E/dV², J/m^6 = N/m^5
#     expt_curv  = curv_ * (10**-60) * (6.242*10**18) # eV / A^6
#
#
#     e_nonx  = energies - contribs @ coefs # len-5 vector
#
#     # also because we only need "a", we only dot the last row with the coef (col) vec
#
#     vander = vstack((ones(len(volumes)),volumes,volumes**2)).T  # 5 x 3
#     vinv   = inv(vander.T @ vander)                             # 3 x 3
#     solver =  vinv @ vander.T                                   # 3 x 5
#     vecs    = solver @ contribs                                        # 3 x 64*5
#     constvec  = solver @ e_nonx                                 # 1 x 3
#
#     curv_vec   = vecs[2]
#     curv_const = expt_curv - constvec[2]
#     lat_vec    = vecs[1] / 2*expt_curv
#     lat_const  = expt_volume -  constvec[1] / 2*expt_curv
#
#     return  curv_vec.tolist(), curv_const, lat_vec.tolist(), lat_const
    # q3 = '''SELECT expt.n_atoms, species.n_atoms, volumes, energies,
    #                 contribs, expt.composition, atomic_contribs,
    #                 atomic_energies, bulk_contribs, bulk_energy,
    #                 expt_cohesive_energy, expt_bm, expt_volume
    #         FROM expt JOIN species ON species=species_id
    #         WHERE calc = %s
    #               AND name REGEXP BINARY(%s)'''
    #
    # q4 = '''SELECT pw,econv,data,a11,a12,a13,a14,a15,msb
    #         FROM calc JOIN functional ON functional=functional_id
    #         WHERE calc_id = %s'''
    #
    # q5 = '''SELECT const_name,val,kind,s,alpha
    #         FROM const'''
    #
    # q7 = '''SELECT expt.n_atoms,expt.composition,symmetry
    #         FROM expt JOIN species ON species=species_id
    #         WHERE calc = %s AND name REGEXP BINARY(%s)'''

        #
        # # Get calculator information
        # #---------------------------
        # pw_,econv_,fx_,a1_,a2_,a3_,a4_,a5_,msb_ = sqlselect(conn,q4,[calc])[0]
        # pw,econv,a11,a12,a13,a14,a15,msb = map(float,[pw_,econv_,a1_,a2_,a3_,a4_,a5_,msb_])
        # fx = loads(fx_)
        # a1s = [a11,a12,a13,a14,a15]
        # # Get input expt information
        # #---------------------------
        #
        # expts = []
        # for n,c,s in sqlselect(conn,q7,[calc,dc]):
        #     expts.append([n,c,s]) # identifying info for an expt, given a calc

        # cons,cnames = [],[]
        #
        # for name,val,kind,s_,alpha in sqlselect(conn,q5):
        #     if name in conweight:
        #         cnames.append(name)
        #         w,s,a,v   = [float(x) if x is not None else None for x in
        #                         [conweight[name],s_,alpha,val]]
        #         cons.append(dict(name=name,weight=w,s=s,alpha=a,kind=kind,val=v))


            # # Get fitting data
            # #-----------------
            # sli = slice(64*i,64*(i+1))
            #
            # fxi = fx[sli]
            #
            # ceX,bmX,lX,ceY,bmY,lY = [], [],[],[],[],[]
            #
            # for n,nsp,vols,engs,xcs,comp,ac,ae,bc,be,ce,bm,vol in sqlselect(conn,q3,[calc,dc]):
            #     a_contribs = {k:array(v)[sli] for k,v in literal_eval(ac).items()}
            #     b_contribs = array(loads(bc)[sli])
            #     contribvec = [xc[sli] for xc in loads(xcs)]
            #     new_cex,new_cey = cohesive(comp,a_contribs,ae,fx,n//nsp,ce,be,b_contribs)
            #     ceX.append(new_cex);ceY.append(new_cey)
            #     new_bmx,new_bmy,new_lx,new_ly = bm_lat(engs,vols,contribvec,float(bm),float(vol),fx)
            #     lX.append(new_lx);bmX.append(new_bmx)
            #     bmY.append(new_ly);lY.append(new_bmy)

            # def write(fi:str,x:Any)->None:
            #     with open(join(root,fi)+'.json','w') as file: dump(x,file)

            # # Write to directory
            # #--------------------
            # write('metadata', md)
            # write('constraints', cons)
            # write('data', [ceX,bmX,lX,ceY,bmY,lY])
            # copyfile(rf,join(root,'runfit.py'))
            # copyfile(sf,join(root,'subfit.sh'))
            # system('chmod 755 '+join(root,'subfit.sh'))
