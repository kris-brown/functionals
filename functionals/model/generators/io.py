# External
from typing import Set as S, Optional as O, Tuple as T
from os.path import join

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Query, Const, GROUP_CONCAT, LEFT,
                   Literal)

from functionals.scripts.io.parse_job import parse_job
from functionals.scripts.io.pop_surf import pop_surf
from functionals.scripts.io.all_mats import all_mats
from functionals.scripts.fit.pop_cons import pop_cons
from functionals.scripts.fit.pop_fitp import pop_fitp

###################################################################
###################################################################
###################################################################
root = '/Users/ksb/scp_tmp/vauto/'
csvpth = '/'+join(*__file__.split('/')[:-4], 'data/%s.csv')
pj_cols = ['stordir', 'runtime', 'pw', 'econv', 'fx', 'contribs', 'energy',
           'int_occupation', 'mag']

mats = elems = set()  # type: S[str]


def io(mod: Model) -> None:

    # Extract tables
    tabs = ['Atoms', 'Functional', "Calc", "Bulkjob", "Bulks", "Fit",
            "Fitparams", "Const", "Surf"]
    Atoms, Fx, Calc, Bulkjob, Bulks, Fit, Fitparams, Constr, Surf = \
        map(mod.get, tabs)

    ###################################################################
    ###################################################################
    ###################################################################

    papb = PyBlock(parse_job, args=[Const(root+'atoms'), Const(elems)],
                   outnames=pj_cols)

    ia = Atoms(insert=True, **{x: papb[x] for x in pj_cols})
    pop_atoms =                               \
        Gen(name='pop_atoms',
            desc='populate Atoms', funcs=[papb], actions=[ia])
    #########################################################################

    pbpb = PyBlock(parse_job, args=[Const(root+'bulks'), Const(mats)],
                   outnames=pj_cols)

    pop_bulkjobs =                               \
        Gen(name='pop_bulkjobs',
            desc='populates BulkJobs',
            funcs=[pbpb],
            actions=[Bulkjob(insert=True, **{x: pbpb[x] for x in pj_cols})])

    #########################################################################

    #########################################################################
    xcols = ['expt_ce', 'expt_bm', 'expt_lat', 'expt_mag']
    ptq = Query(dict(b=Bulks.id(), n=Bulks['name']()))

    def parse_csv(root: str, mat: str
                  ) -> T[O[float], O[float], O[float], O[float]]:
        import csv
        if mat in ['PdC', "PdN", "ScC", "NbC", "NbN", 'Re', 'Co', 'Os', 'Tc', 'Tl', 'Y', 'Zr']:
            return None, None, None, None
        bad = ['Re', 'Os', 'Ir', 'Hf', 'Pb', 'W', 'Zn', 'Ta', 'Pd', 'Pt']
        with open(root % 'expt', 'r') as f:
            r = csv.reader(f)
            for ro in r:
                if ro[0] == mat:
                    x, b, z, m = [float(x) if x else None for x in ro[1:]]
                    if any(mm in mat for mm in bad):
                        x = None
                    return x, b, z, m
            print('\n\n\n', mat)
            import pdb
            pdb.set_trace()
            raise ValueError(mat)

    ptpb = PyBlock(parse_csv,
                   args=[Const(csvpth), ptq['n']],
                   outnames=xcols)

    pop_expt = Gen(name='pop_expt',
                   desc='tran https://aip.scitation.org/doi/10.1063/1.4948636',
                   query=ptq, funcs=[ptpb],
                   actions=[Bulks(bulks=ptq['b'],
                                  **{x: ptpb[x] for x in xcols})])
    ##########################################################################

    pccols = ['name', 'abdata', 'kind', 'points', 'val', 'func']
    pcq = Query(dict(x=GROUP_CONCAT(Fitparams['consts'](), delim='$')),
                aggcols=[LEFT(Fitparams['consts'](), Literal(1))])
    pcpb = PyBlock(pop_cons, args=[pcq['x']], outnames=pccols)
    popc = \
        Gen(name='popcon', funcs=[pcpb], tags=['fit'], query=pcq,
            actions=[Constr(insert=True, **{x: pcpb[x] for x in pccols})])
    ##########################################################################

    pfcols = ['ce_scale', 'bm_scale', 'lc_scale', 'consts']
    pfpb = PyBlock(pop_fitp, outnames=pfcols)

    popfp = \
        Gen(name='popfp', tags=['fit'],
            funcs=[pfpb], actions=[Fitparams(insert=True,
                                             **{x: pfpb[x] for x in pfcols})])

    amq = Query(dict(c=Calc.id()))
    ampb = PyBlock(all_mats, args=[Const(csvpth)])
    allmat =    \
        Gen(name='allmat',
            desc='Goes through CSV files to get list of all mats with data',
            funcs=[ampb], query=amq,
            actions=[Calc(calc=amq['c'], allmat=ampb['out'])])
    #######################################################################
    surfcol = ['mat', 'xc', 'ontop', 'hollow']
    spb = PyBlock(pop_surf, outnames=surfcol)
    surf = \
        Gen(name='pop_surf',
            funcs=[spb], tags=['surf'],
            actions=[Surf(insert=True, **{x: spb[x] for x in surfcol})])
    #######################################################################
    #######################################################################
    #######################################################################

    gens = [pop_atoms, pop_bulkjobs, pop_expt, popc, popfp, allmat, surf]

    mod.add(gens)
