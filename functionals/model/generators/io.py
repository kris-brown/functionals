# External
from typing import Optional as O, Tuple as T
from os.path import join

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Query, Const)

from functionals.scripts.io.parse_job import parse_job
from functionals.scripts.io.parse_bulks import parse_bulks
from functionals.scripts.io.pop_surf import pop_surf
from functionals.scripts.io.all_mats import all_mats
from functionals.scripts.io.scheffdiff import scheffdiff
from functionals.scripts.fit.pop_fitp import pop_fitp

###################################################################
###################################################################
###################################################################
root = '/Users/ksb/scp_tmp/vauto/'
csvpth = '/' + join(*__file__.split('/')[:-4], 'data/%s.csv')


def io(mod: Model) -> None:

    # Extract tables
    tabs = ['Atoms', "Calc", "Bulks", "Fit",
            "Fitparams", "Surf"]
    Atoms, Calc, Bulks, Fit, Fitparams, Surf = \
        map(mod.get, tabs)

    ###################################################################
    ###################################################################
    ###################################################################
    pj_cols = ['stordir', 'contribs', 'energy', 'mag', 'err']

    papb = PyBlock(parse_job, args=[Const(root + 'atoms')],
                   outnames=pj_cols)
    papb2 = PyBlock(lambda ns: [n.split('/')[-2]
                                for n in ns], args=[papb['stordir']])
    iac = Calc(insert=True, name=papb2['out'])
    ia = Atoms(insert=True, calc=iac, calcname=papb2['out'],
               **{x: papb[x] for x in pj_cols if x != 'mag'})
    pop_atoms =                               \
        Gen(name='pop_atoms',
            desc='populate Atoms', transforms=[papb, papb2], loads=[ia])
    #########################################################################
    pb_cols = ['stordir', 'name', 'calcname', 'eng', 'contrib',
               'mag', 'volumes', 'energies', 'contribs', 'eos', 'err']

    pbpb = PyBlock(parse_bulks, args=[Const(root + 'bulks')], outnames=pb_cols)
    ibc = Calc(insert=True, name=pbpb['calcname'])
    ib = Bulks(insert=True, calc=ibc, **{x: pbpb[x] for x in pb_cols})

    pop_bulkjobs = \
        Gen(name='pop_bulkjobs',
            desc='populates Bulks',
            transforms=[pbpb], loads=[ib])

    #########################################################################
    xcols = ['expt_ce', 'expt_bm', 'expt_lat', 'expt_mag']
    ptq = Query(dict(b=Bulks.id(), n=Bulks['name']()))

    def pop_exp(root: str, mat: str
                ) -> T[O[float], O[float], O[float], O[float]]:
        import csv
        with open(root % 'expt', 'r') as f:
            r = csv.reader(f)
            for ro in r:
                if ro[0] == mat:
                    x, b, z, m = [float(x) if x else None for x in ro[1:]]
                    return x, b, z, m
        raise ValueError("No Expt data!")

    ptpb = PyBlock(pop_exp,
                   args=[Const(csvpth), ptq['n']],
                   outnames=xcols)

    pop_expt = Gen(name='pop_expt',
                   desc='tran https://aip.scitation.org/doi/10.1063/1.4948636',
                   query=ptq, transforms=[ptpb],
                   loads=[Bulks(bulks=ptq['b'],
                                **{x: ptpb[x] for x in xcols})])
    ##########################################################################

    # pccols = ['name', 'abdata', 'kind', 'points', 'val', 'func']
    # pcq = Query(dict(x=GROUP_CONCAT(Fitparams['consts'](), delim='$')),
    #             aggcols=[LEFT(Fitparams['consts'](), Literal(1))])
    # pcpb = PyBlock(pop_cons, args=[pcq['x']], outnames=pccols)
    # popc = \
    #     Gen(name='popcon', transforms=[pcpb], tags=['fit'], query=pcq,
    #         loads=[Constr(insert=True, **{x: pcpb[x] for x in pccols})])
    ##########################################################################

    pfcols = ['ce_scale', 'bm_scale', 'lc_scale']  # , 'consts', 'abdata']
    pfpb = PyBlock(pop_fitp, outnames=pfcols)

    popfp = \
        Gen(name='popfp', tags=['fit'],
            transforms=[pfpb], loads=[Fitparams(insert=True,
                                                **{x: pfpb[x] for x in pfcols})])

    # amq = Query(dict(c=Calc.id()))
    # ampb = PyBlock(all_mats, args=[Const(csvpth)])
    # allmat =    \
    #     Gen(name='allmat',
    #         desc='Goes through CSV files to get list of all mats with data',
    #         transforms=[ampb], query=amq,
    #         loads=[Calc(calc=amq['c'], allmat=ampb['out'])])

    #######################################################################
    surfcol = ['mat', 'xc', 'ontop', 'hollow']
    spb = PyBlock(pop_surf, outnames=surfcol)
    surf = \
        Gen(name='pop_surf',
            transforms=[spb], tags=['surf'],
            loads=[Surf(insert=True, **{x: spb[x] for x in surfcol})])
    #######################################################################
    bq = Query(dict(c=Calc.id(), n=Calc['name']()))

    def bf(n: str) -> O[str]:
        import os
        root = '/Users/ksb/functionals/data/beefs'
        for b in os.listdir(root):
            if n == b[:-5]:
                with open(os.path.join(root, b), 'r') as f:
                    return f.read()
        return None
    bpb = PyBlock(bf, args=[bq['n']])
    beef = Gen(name='beef', query=bq, transforms=[bpb],
               loads=[Calc(calc=bq['c'], data=bpb['out'])])
    #######################################################################
    scheffout = ['dscheff_' + k for k in ['ce', 'bm', 'lat']]
    errs = [Bulks['err_' + k]() for k in ['ce', 'bm', 'lat']]
    schq = Query(dict(z=Bulks.id(), n=Bulks['name'](), c=Bulks['calcname'](),
                      e=errs[0], b=errs[1], t=errs[2]),
                 opt_attr=errs)

    schpb = PyBlock(scheffdiff, args=[Const(csvpth % "scheff"),
                                      *[schq[k] for k in 'ncebt']],
                    outnames=scheffout)
    scheff = Gen(name='scheffdiff', query=schq, transforms=[schpb],
                 loads=[Bulks(bulks=schq['z'],
                              **{x: schpb[x] for x in scheffout})])
    #######################################################################
    #######################################################################

    gens = [pop_atoms, pop_bulkjobs, pop_expt, popfp, surf, beef,
            scheff]

    mod.add(gens)
