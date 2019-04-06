# External
from os import environ
from os.path import join
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
    tabs = ['Fitparams','Const','Atoms','Functional',"Job","Calc","Bulks",'Fit']
    Fitparams,Cnst, Atoms,Fx, Job, Calc, Bulks, Fit = map(mod.get,tabs)

    ################################################################################
    ################################################################################
    ################################################################################

    fitcols = ['constden','consts','reg', 'bm_weight','lat_weight']

    parampth = root('data/fitparams.csv')

    pfpb = PyBlock(parse_fitparams,
                   args     = [Const(parampth)],
                   outnames = fitcols)

    pop_fitparams =                                                             \
        Gen(name    = 'pop_fitparams',
            desc    = 'Looks in FITPATH for fitting specifications',
            funcs   = [pfpb],
            tags    = ['fit'],
            actions = [Fitparams(insert = True,
                           **{x:pfpb[x] for x in fitcols})])

    ################################################################################

    constpth = root('data/constraints.csv')
    constcols = ['const_name','description','s','alpha','val','kind']

    pcpb      = PyBlock(parse_const, args = [Const(constpth)], outnames = constcols)

    pop_constraint =                                                            \
        Gen(name    = 'pop_constraint',
            desc    = 'populate constraints',
            funcs   = [pcpb],
            tags    = ['fit'],
            actions = [Cnst(insert = True, **{x:pcpb[x] for x in constcols})])


    ############################################################################

    acols  = ['contribs','energy','int_occupation','num','true_mag','mag']
    jcols  = ['stordir','name','runtime','pw','econv','fx'] # job cols
    pacols = jcols + acols
    apth  = '/Users/ksb/scp_tmp/vauto/atoms'
    papb  = PyBlock(parse_atoms, args = [Const(apth)], outnames = pacols)

    ifun  = Fx(insert = True, data   = papb['fx'])

    ijob  = Job(insert  = True,
                calc    = Calc(insert     = True,
                               pw         = papb['pw'],
                               econv      = papb['econv'],
                               functional = ifun),
                stordir = papb['stordir'])

    pop_atoms =                                                                 \
        Gen(name    = 'pop_atoms',
            desc    = 'populate Atoms',
            funcs   = [papb],
            tags    = ['io'],
            actions = [Atoms(insert=True,
                             job   = ijob,
                             name  = papb['name'],
                             **{x:papb[x] for x in acols})])

    ############################################################################

    bcols  = ['contribs','energies','volumes','n_atoms','n_elems','composition','img',
              'eosbm', 'lattice','completed','success']
    bpth   = '/Users/ksb/scp_tmp/vauto/bulks'
    pbcols = jcols + bcols
    pbpb   = PyBlock(parse_bulks, args = [Const(bpth)], outnames = pbcols)

    ifun_  = Fx(insert = True, data   = pbpb['fx'])

    ijob_  = Job(insert  = True,
                 calc    = Calc(insert     = True,
                               pw         = pbpb['pw'],
                               econv      = pbpb['econv'],
                               functional = ifun_),
                 stordir = pbpb['stordir'])

    pop_bulks =                                                                 \
        Gen(name = 'pop_bulks',
            desc = 'populates Bulks',
            funcs = [pbpb],
            tags  = ['io'],
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
    #a1msb = ['a11','a12','a13','a14','a15','msb']

    fcols = ['name','a1','pth','timestamp','decay','nsteps','steps','result','runtime']

    fitpth   = join(environ['HOME'],'scp_tmp/fitresult')

    fpb = PyBlock(parse_fit,
                  args     = [Const(fitpth)],
                  outnames = fcols + ['pw','econv','data'] + fitcols)

    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'looks for completed fitting jobs - grabs identifying info',
            funcs   = [fpb],
            tags    = ['io','fit'],
            actions = [Fit(insert = True,
                           **{x:fpb[x] for x in fcols},
                           calc   = Calc(insert = True,
                                        **{x:fpb[x] for x in ['pw','econv']},
                                        functional = Fx(insert = True,
                                                        data   = fpb['data'])),
                           fitparams = Fitparams(insert=True,
                                                 **{x:fpb[x] for x in fitcols}))])



    ############################################################################
    ############################################################################
    ############################################################################

    gens = [pop_fitparams,pop_constraint,pop_atoms,pop_bulks,pop_expt,pop_fit]

    mod.add(gens)
