from os.path import join

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, Expr, Literal as Lit, COALESCE,
                   Text, JPath, GROUP_CONCAT, MAX, AND, CONVERT, SUBSELECT)

from functionals.scripts.fit.db_data import db_data
from functionals.scripts.fit.a_ce import a_ce
from functionals.scripts.fit.fit_vectors import fit_vectors
from functionals.scripts.fit.runfit import runfit

############################################################################
############################################################################
############################################################################
functionals_db = '/'+join(*__file__.split('/')[:-3], 'data/functionals.json')


def fit(mod: Model) -> None:
    # Extract tables

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

    refpth = JPath("refs", [refs__bulks])
    beefpth = JPath('calc', [bulks__calc])
    aceq = Query(exprs=dict(
        b=Bulks.id(),
        r=GROUP_CONCAT(Refs['energy'](refpth)),
        c=GROUP_CONCAT(Refs['contribs'](refpth), delim='$'),
        e=Bulks['eng'](),
        x=Bulks['contribs'](),
        n=Bulks['n_atoms'](),
        f=MAX(Calc['data'](beefpth)),
        m=Bulks['name']()),
        aggcols=[Bulks.id()],
        basis=[Bulks])

    acepb = PyBlock(a_ce, args=[aceq[x] for x in 'rcexnfm'],)

    ace =                                                          \
        Gen(name='ace',
            desc='Computes A matrix vectors to yield cohesive energy',
            query=aceq,
            funcs=[acepb], tags=['fit'],
            actions=[Bulks(bulks=aceq['b'], ab_ce=acepb['out'])])

    ########################################################################
    vecout = ['ab_'+x for x in ['bm', 'vol']]

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
                            n=Bulks['name']()),
                 basis=[Bulks],
                 opt_attr=[Bulks['expt_bm'](), Bulks['expt_vol']()])

    vecpb = PyBlock(fit_vectors,
                    args=[vecq[x] for x in 'evcfon'],
                    outnames=vecout)

    vecs =                                                            \
        Gen(name='vecs',
            desc='Computes Ax+b=y matrices/vectors',
            query=vecq, funcs=[vecpb],
            tags=['fit'],
            actions=[Bulks(bulks=vecq['b'], **{x: vecpb[x] for x in vecout})])

    ########################################################################

    bpth = JPath('bulks', [bulks__calc])
    names = ['ab_ce', 'ab_bm', 'ab_vol', 'expt_ce', 'expt_bm', 'expt_vol',
             'name', 'ce']
    fdgbs = dict(zip('abcdefgh', names))

    def null_to_empty(x: Expr) -> Expr:
        '''Cast attribute to string - if null, return empty string.'''
        return COALESCE([CONVERT(x, Text()), Lit('')])

    fddict = {k: GROUP_CONCAT(null_to_empty(Bulks[v](bpth)), delim='$')
              for k, v in fdgbs.items()}

    fdq = Query(exprs=dict(z=Calc.id(), **fddict),
                basis=[Calc],
                aggcols=[Calc.id()],
                opt_attr=[Bulks[v](bpth) for v in fdgbs.values()],
                constr=AND([Calc['beef']()]))

    fdpb = PyBlock(db_data, args=[fdq[x] for x in 'abcdefgh'])

    fitdata = Gen(
        name='fitdata',
        desc='Assembles all BEEF data into objects for fitting',
        query=fdq,
        funcs=[fdpb],
        tags=['fit'],
        actions=[Calc(calc=fdq['z'], fitdata=fdpb['out'])]
    )

    ########################################################################

    funmetrics = ['mae_%s' % (x) for x in ['ce', 'bm', 'lat']]

    gc = SUBSELECT("""CONCAT('{',
string_agg(CONCAT('''', name, ''':[', abdata, ',''', kind, ''','
                  , points, ']'), ','),
'}')""", 'const', 'true')
    rfq = Query(dict(c=Calc.id(), d=Calc['fitdata'](), n=gc,
                     o=Fitparams['consts'](),
                     p=Fitparams.id(), e=Fitparams['ce_scale'](),
                     b=Fitparams['bm_scale'](), v=Fitparams['lc_scale']()),
                basis=[Fitparams, Calc])

    rfout = ['x', 'cv'] + funmetrics
    rfpb = PyBlock(runfit, args=[rfq[x] for x in 'debvon'], outnames=rfout)
    rfit = Gen(name='runfit',
               desc='', tags=['fit'],  # 'parallel'],
               query=rfq, funcs=[rfpb],
               actions=[Fit(insert=True, calc=rfq['c'], fitparams=rfq['p'],
                            **{x: rfpb[x] for x in rfout})])
    ########################################################################
    ########################################################################
    ########################################################################
    gens = [ace, vecs, fitdata, rfit]
    mod.add(gens)
