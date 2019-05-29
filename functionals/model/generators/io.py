# External
from typing import Tuple as T
from os import environ
from os.path import join
from json import load,dumps

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Env, Query, Const, Import, defaultEnv,
                    EQ, Literal as Lit, Arg, JPath, CONCAT, NULL, LIKE, NOT, AND)

from functionals.scripts.io.parse_fitparams        import parse_fitparams
from functionals.scripts.io.parse_const            import parse_const
from functionals.scripts.io.parse_atoms            import parse_atoms
from functionals.scripts.io.parse_bulks            import parse_bulks
from functionals.scripts.io.parse_keld             import parse_keld
from functionals.scripts.io.parse_fit              import parse_fit

################################################################################
################################################################################
################################################################################
def root(x:str)->str: return join(environ['FUNCTIONALS_ROOT'],x)

def io(mod : Model) -> None:

    # Extract tables
    tabs = ['Fitparams','Atoms','Functional',"Job","Calc","Bulks",'Fit']
    Fitparams, Atoms,Fx, Job, Calc, Bulks, Fit = map(mod.get,tabs)

    ################################################################################
    ################################################################################
    ################################################################################

    fitcols = ['consts','reg', 'ce_scale','bm_scale','lc_scale']

    parampth = root('data/fitparams.csv')

    pfpb = PyBlock(parse_fitparams,
                   args     = [Const(parampth)],
                   outnames = fitcols)

    pop_fitparams =                                                             \
        Gen(name    = 'pop_fitparams',
            desc    = 'Looks in FITPATH for fitting specifications',
            funcs   = [pfpb],
            tags    = ['fit','io'],
            actions = [Fitparams(insert = True,
                           **{x:pfpb[x] for x in fitcols})])


    ############################################################################

    acols  = ['contribs','energy','int_occupation','num','true_mag','mag']
    jcols  = ['stordir','name','runtime','pw','fx'] # job cols
    pacols = jcols + acols
    apth  = '/Users/ksb/scp_tmp/vauto/atoms'
    papb  = PyBlock(parse_atoms, args = [Const(apth)], outnames = pacols)

    ifun  = Fx(insert = True, data   = papb['fx'])

    ijob  = Job(insert  = True,
                runtime = papb['runtime'],
                calc    = Calc(insert     = True,
                               pw         = papb['pw'],
                               #econv      = papb['econv'],
                               functional = ifun),
                stordir = papb['stordir'])

    pop_atoms =                                                                 \
        Gen(name    = 'pop_atoms',
            desc    = 'populate Atoms',
            funcs   = [papb],
            tags    = [],#['io'],
            actions = [Atoms(insert=True,
                             job   = ijob,
                             name  = papb['name'],
                             **{x:papb[x] for x in acols})])

    ############################################################################

    bcols  = ['alleng','allvol','allcontrib','n_atoms','n_elems','composition','img',
              'eosbm', 'lattice','strain_low','strain_hi','morejobs','incomplete','success']
    bpth   = '/Users/ksb/scp_tmp/vauto/bulks'
    pbcols = jcols + bcols
    pbpb   = PyBlock(parse_bulks, args = [Const(bpth)], outnames = pbcols)

    ifun_  = Fx(insert = True, data   = pbpb['fx'])

    ijob_  = Job(insert  = True,
                 runtime    = pbpb['runtime'],
                 calc    = Calc(insert     = True,
                                pw         = pbpb['pw'],
                                functional = ifun_),
                 stordir = pbpb['stordir'])

    pop_bulks =                                                                 \
        Gen(name = 'pop_bulks',
            desc = 'populates Bulks',
            funcs = [pbpb],
            tags  = [],#['io'],
            actions = [Bulks(insert=True,job=ijob_,name=pbpb['name'],
                            **{x:pbpb[x] for x in bcols})])

    ############################################################################
    xcols    = ['expt_ce', 'expt_bm', 'expt_l', 'expt_mag']
    keldpath = environ['KELDPATH']

    peq = Query(dict(b = Bulks.id(), n = Bulks['name']()))

    pepb = PyBlock(parse_keld,
                   args     = [peq['n'],Const(keldpath)],
                   outnames = xcols)

    pop_expt =                                                                  \
        Gen(name    = 'pop_expt',
            desc    = 'Populates experimental data cols in Bulks',
            query   = peq,
            funcs   = [pepb],
            actions = [Bulks(bulks=peq['b'],**{x:pepb[x] for x in xcols})])

    ############################################################################

    fcols = ['name','pth']

    fitpth   = join(environ['FUNCTIONALS_ROOT'],'data/fit')

    fpb = PyBlock(parse_fit,
                  args     = [Const(fitpth)],
                  outnames = fcols + ['pw','fdata'] + fitcols)

    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'looks for completed fitting jobs - grabs identifying info',
            funcs   = [fpb],
            tags    = ['io','fit'],
            actions = [Fit(insert = True,
                           **{x:fpb[x] for x in fcols},
                           calc   = Calc(insert = True,
                                         pw = fpb['pw'],
                                        functional = Fx(insert = True,
                                                        data   = fpb['fdata'])),
                           fitparams = Fitparams(insert=True,
                                                 **{x:fpb[x] for x in fitcols}))])



    ############################################################################
    ptcols = ['traj%d'%i for i in range(5)]

    ptq = Query(dict(f=Fit.id(),p=Fit['pth']()))

    def ptfun(pth:str)->T[str,str,str,str,str]:
        with open(join(pth,'result.json'),'r') as f: a,b,c,d,e = map(dumps,load(f))
        return a,b,c,d,e

    ptpb = PyBlock(ptfun,env=defaultEnv+Env(Import('os.path','join')),args=[ptq['p']],outnames=ptcols)

    pop_traj = \
        Gen(name = 'pop_traj',
            desc = 'loads trajectories of fits',
            query = ptq,
            funcs = [ptpb],
            actions = [Fit(fit=ptq['f'],**{x:ptpb[x] for x in ptcols})])
    ############################################################################
    ############################################################################
    ############################################################################

    gens = [pop_fitparams,pop_atoms,pop_bulks,pop_expt,pop_fit,pop_traj]

    mod.add(gens)
