# External
from typing import Tuple as T, Set as S
from os import environ
from os.path import join
from csv import reader
from json import load,dumps

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Env, Query, Const, Import, defaultEnv,
                    EQ, Literal as Lit, Arg, JPath, CONCAT, NULL, LIKE, NOT,
                    AND, GROUP_CONCAT, MAX)

from functionals.scripts.io.parse_const          import parse_const
from functionals.scripts.io.parse_csv            import parse_csv
from functionals.scripts.io.parse_job            import parse_job
from functionals.scripts.io.parse_fit            import parse_fit

################################################################################
################################################################################
################################################################################
root = '/Users/ksb/scp_tmp/vauto/'
csvpth = '/'+join(*__file__.split('/')[:-4], 'data/%s.csv')
pj_cols = ['stordir','runtime','pw','econv','fx','contribs','energy', 'int_occupation', 'mag']

# with open(csvpth,'r') as f:
#     csvr = reader(f)
#     mats = set(['/%s_'%row[0] for row in csvr if 'hcp' != row[1]])
# elems = {'/%s$'%x for x in {'Ag','Al','As','N','P','Au','B','Ba','Be','C','Ca','O','Cd','Co','Cu',
#          'Fe','Ge','As','Ga','Ge','In','Ir','K','Li','Cl','F','H','Mg','S','Mo',
#          'Na','Nb','Ni','Os','Pb','Pd','Pt','Rb','Rh','Ru','Sc','Si','Sn','Sr',
#          'Ta','Ti','V','W','Zn','Zr'}}
mats = elems = set() # type: S[str]


def io(mod : Model) -> None:

    # Extract tables
    tabs = ['Atoms','Functional',"Calc","Bulkjob", "Bulks", "Fit", "Fitparams"]
    Atoms, Fx, Calc, Bulkjob, Bulks, Fit, Fitparams = map(mod.get,tabs)

    ################################################################################
    ################################################################################
    ################################################################################

    papb  = PyBlock(parse_job, args = [Const(root+'atoms'),Const(elems)], outnames = pj_cols)

    ia = Atoms(insert=True, **{x:papb[x] for x in pj_cols})
    pop_atoms =                                                                 \
        Gen(name    = 'pop_atoms',
            desc    = 'populate Atoms',
            funcs   = [papb],
            actions = [Atoms(insert=True, **{x:papb[x] for x in pj_cols})])
    ############################################################################

    pbpb   = PyBlock(parse_job, args = [Const(root+'bulks'),Const(mats)], outnames = pj_cols)


    pop_bulkjobs =                                                                 \
        Gen(name = 'pop_bulkjobs',
            desc = 'populates BulkJobs',
            funcs = [pbpb],
            actions = [Bulkjob(insert=True, **{x:pbpb[x] for x in pj_cols})])

    ############################################################################

    ############################################################################
    xcols = ['expt_ce', 'expt_bm', 'expt_l','expt_vol', 'expt_mag']
    ptq = Query(dict(b = Bulks.id(), n = Bulks['name'](), r=Bulks['volrat']()))

    ptpb = PyBlock(parse_csv,
                   args = [Const(csvpth), ptq['n'], ptq['r']],
                   outnames = xcols)

    pop_expt = Gen(name = 'pop_expt',
                   desc='tran https://aip.scitation.org/doi/10.1063/1.4948636',
                   query= ptq, funcs=[ptpb],
                   actions = [Bulks(bulks=ptq['b'],**{x:ptpb[x] for x in xcols})])
    ############################################################################

    fcols = ['name','pth','done']
    fitcols = ['consts','reg', 'ce_scale','bm_scale','lc_scale','mag_scale']

    env = Env([Import('os.path',['join','getsize']),
               Import('os',['listdir'])]) + defaultEnv
    fpb = PyBlock(parse_fit,
                  env = env,
                  args     = [Const(root+'fit')],
                  outnames = fcols + ['pw','fdata'] + fitcols)

    icalc = Calc(insert=True,
                 pw=fpb['pw'],
                 functional=Fx(insert=True, data=fpb['fdata']))

    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'looks for completed fitting jobs - grabs identifying info',
            funcs   = [fpb],
            tags    = ['io','fit'],
            actions = [Fit(insert = True,calc=icalc,
                           **{x:fpb[x] for x in fcols},
                           fitparams = Fitparams(insert=True,
                                                 **{x:fpb[x] for x in fitcols}))])



    ############################################################################
    ############################################################################
    ############################################################################

    gens = [pop_atoms, pop_bulkjobs, pop_fit, pop_expt]

    mod.add(gens)
