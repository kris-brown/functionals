from typing import List as L, Tuple as T

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, COALESCE,
                   Text, GROUP_CONCAT, CONCAT, MAX, AND, CONVERT, GT, LEFT,
                   ABS, EQ, Varchar, Const)

from functionals.scripts.fit.db_data import db_data
from functionals.scripts.fit.a_ce import a_ce
from functionals.scripts.fit.fit_vectors import fit_vectors
from functionals.scripts.fit.runfit import runfit

############################################################################
############################################################################
############################################################################


def fit(mod: Model) -> None:

    # Extract tables
    tabs = ['fit', 'calc', 'fitparams', 'bulks', 'refs']
    Fit, Calc, Fitparams, Bulks, Refs = map(mod.get, tabs)

    # Extract rels
    fit__calc, bulks__calc, refs__bulks, fit__fitparams = \
        map(mod.get_rel, [Fit.r('calc'),
                          Bulks.r("calc"), Refs.r('bulks'),
                          Fit.r('fitparams')])
    ########################################################################
    ########################################################################

    refpth = mod.make_path("refs", [refs__bulks])
    beefpth = mod.make_path('calc', [bulks__calc])
    aceq = Query(exprs=dict(
        b=Bulks.id(),
        r=GROUP_CONCAT(Refs['energy'](refpth)),
        c=GROUP_CONCAT(Refs['contribs'](refpth), delim='$'),
        e=Bulks['eng'](),
        x=Bulks['contribs'](),
        n=Bulks['n_atoms'](),
        f=MAX(Calc['data'](beefpth)),
        m=Bulks['name'](),
        t=Bulks['ce']()),
        aggcols=[Bulks.id()],
        basis=[Bulks])

    acepb = PyBlock(a_ce, args=[aceq[x] for x in 'rcexnfmt'],)

    ace =                                                          \
        Gen(name='ace',
            desc='Computes A matrix vectors to yield cohesive energy',
            query=aceq,
            transforms=[acepb], tags=['fit'],
            loads=[Bulks(bulks=aceq['b'], ab_ce=acepb['out'])])

    ########################################################################
    from functionals.fit.functional import FromMatrix
    import json
    # ms2c = Const(json.dumps(FromMatrix.frompath('ms2').x.tolist()))
    ms2c = Const([json.dumps(FromMatrix.frompath('ms2').x.tolist()),
                  json.dumps(FromMatrix.frompath('pbesol').x.tolist()),
                  json.dumps(FromMatrix.frompath('7500').x.tolist()),
                  ])

    c0 = Const([0, 0, 0])

    ms2q = Query(dict(c=Calc.id()), constr=EQ(Calc['name'](), Lit('ms2')))
    fact = Fitparams(insert=True, ce_scale=c0, bm_scale=c0, lc_scale=c0)
    ms2act = Fit(x=ms2c, step=c0, fitparams=fact, calc=ms2q['c'], insert=True)

    ms2fit =                                                          \
        Gen(name='ms2fit',
            desc='Makes a fake "fit results" equal to MS2/PBE/7500',
            query=ms2q, tags=['fit'], loads=[ms2act])

    ########################################################################
    vecout = ['ab_' + x for x in ['bm', 'vol']]

    aceq = Query(exprs=dict(b=Bulks.id(),
                            r=GROUP_CONCAT(Refs['energy'](refpth)),
                            c=GROUP_CONCAT(Refs['contribs'](refpth),
                                           delim='$'),
                            x=Bulks['contribs'](),
                            e=Bulks['energies'](),
                            n=Bulks['n_atoms'](),
                            f=MAX(Calc['data'](beefpth))),
                 aggcols=[Bulks.id()],
                 basis=[Bulks])

    vecq = Query(exprs=dict(b=Bulks.id(),
                            e=Bulks['energies'](),
                            v=Bulks['volumes'](),
                            c=Bulks['contribs'](),
                            f=Calc['data'](beefpth),
                            o=Bulks['vol'](),
                            n=Bulks['name'](),
                            m=Bulks['bm'](),
                            p=Bulks['prim'](),
                            r=Bulks['volprimrat']()),
                 basis=[Bulks],
                 opt_attr=[Bulks['expt_bm'](), Bulks['expt_vol']()],
                 # constr=NOT(Bulks['irregular']())
                 )

    vecpb = PyBlock(fit_vectors,
                    args=[vecq[x] for x in 'evcfonmrp'],
                    outnames=vecout)

    vecs =                                                            \
        Gen(name='vecs',
            desc='Computes Ax+b=y matrices/vectors for BM and volume',
            query=vecq, transforms=[vecpb], tags=['fit'],
            loads=[Bulks(bulks=vecq['b'], **{x: vecpb[x] for x in vecout})])

    ########################################################################

    ########################################################################

    bpth = mod.make_path('bulks', [bulks__calc])
    names = ['ab_ce', 'ab_bm', 'ab_vol', 'expt_ce', 'expt_bm', 'expt_vol',
             'name', 'ce', 'volrat']
    fdgbs = dict(zip('abcdefghr', names))

    def null_to_empty(x: Expr) -> Expr:
        '''Cast attribute to string - if null, return empty string.'''
        return COALESCE([CONVERT(x, Text()), Lit('')])

    fddict = {k: GROUP_CONCAT(null_to_empty(Bulks[v](bpth)), delim='$')
              for k, v in fdgbs.items()}

    fdq = Query(exprs=dict(z=Calc.id(), n=Calc['name'](), **fddict),
                basis=[Calc],
                aggcols=[Calc.id()],
                opt_attr=[Bulks[v](bpth) for v in fdgbs.values()],
                constr=AND([Calc['beef']()]))

    fdpb = PyBlock(db_data, args=[fdq[x] for x in 'abcdefghrn'])

    fitdata = Gen(
        name='fitdata',
        desc='Assembles all BEEF data into objects for fitting',
        query=fdq,
        transforms=[fdpb],
        tags=['fit'],
        loads=[Calc(calc=fdq['z'], fitdata=fdpb['out'])]
    )

    ########################################################################

    funmetrics = [y + 'mae_%s' % (x)
                  for x in ['ce', 'bm', 'lat', 'vol', 'mag']
                  for y in ['', 'rel']]

    rfq = Query(dict(c=Calc.id(), d=Calc['fitdata'](),
                     o=Calc['data'](),
                     x=LEFT(CONVERT(ABS(Fit.id()), Varchar()), Lit(4)),
                     p=Fitparams.id(), e=Fitparams['ce_scale'](),
                     b=Fitparams['bm_scale'](), v=Fitparams['lc_scale']()),
                basis=[Fitparams, Calc],
                constr=AND([GT(Fitparams['ce_scale'](), Lit(0)),
                            EQ(Calc['n_missing'](), Lit(0))]))

    rfout = ['x', 'step', 'cv'] + funmetrics + ['bump', 'err']
    rfpb = PyBlock(runfit, args=[rfq[x] for x in 'debvo'], outnames=rfout)
    rfit = Gen(name='runfit',
               desc='', tags=['fit'],  # 'parallel'],
               query=rfq, transforms=[rfpb],
               loads=[Fit(insert=True, calc=rfq['c'], code=rfq['x'],
                          fitparams=rfq['p'],
                          **{x: rfpb[x] for x in rfout})])
    ########################################################################

    fpq = Query(dict(d=GROUP_CONCAT(
        CONCAT([Calc['name'](), Lit('$'),
                CONVERT(Calc.id(), Varchar()), Lit('$'), Calc['data']()]),
        delim='|'),
        f=Fitparams.id()),
        aggcols=[Fitparams.id()],
        basis=[Calc, Fitparams], constr=EQ(Fitparams['ce_scale'](), Lit(0)))

    def fpfn(rawdata: str, fpid: int
             ) -> T[L[int], L[str], L[str], L[int], int]:
        fxs = {}
        for fx in rawdata.split('|'):
            a, b, c = fx.split('$')
            fxs[a] = (b, c)
        pairs = [('ms2', 'ms2'), ('7500', 'ms2'), ('5558', '7500'),
                 ('2751', '5558')]
        steps = list(range(len(pairs)))  # identifying attr, so make different
        names = ['{1} -> {0}'.format(y, x) for x, y in pairs]
        xs = [fxs[x][1] for x, _ in pairs]
        calcs = [int(fxs[x][0]) for _, x in pairs]
        return steps, names, xs, calcs, fpid

    fpout = ['step', 'name', 'x', 'calc', 'fitparams']
    fppb = PyBlock(fpfn, args=[fpq['d'], fpq['f']], outnames=fpout)
    fitpairs = Gen(
        name='fitpairs',
        desc='Inserts fits that declare one fx was derived from another',
        query=fpq, tags=['fit'], transforms=[fppb],
        loads=[Fit(insert=True, **{x: fppb[x] for x in fpout})])
    ########################################################################
    ########################################################################
    ########################################################################
    gens = [ace, vecs, fitdata, rfit, ms2fit, fitpairs]
    mod.add(gens)
