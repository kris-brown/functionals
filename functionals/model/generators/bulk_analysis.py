# Internal
from dbgen import (Model, Gen, Query, PyBlock, NOT, AND,
                   Literal as Lit, CASE, OR, LIKE, NULL,
                   GT, ABS)

from functionals.scripts.load.analyze_bulks import analyze_bulks
from functionals.scripts.load.analyze_opt_bulks import analyze_opt_bulks


"""
Generators purely related to analysis of Bulks and Bulkjobs.
"""


def bulk_analysis(mod: Model) -> None:
    tabs = ['Bulks', 'Calc']
    Bulks, Calc = map(mod.get, tabs)
    ########################################################################
    ########################################################################
    ########################################################################

    def get_struct(mat: str) -> str:
        from functionals.CLI.submit import matdata
        return str(matdata[mat].struct) if mat in matdata else ''

    sq = Query(dict(b=Bulks.id(), n=Bulks['name']()))
    sf = PyBlock(get_struct, args=[sq['n']])
    struct = \
        Gen(name='struct', query=sq, funcs=[sf],
            actions=[Bulks(bulks=sq['b'], struct=sf['out'])])
    ########################################################################
    ########################################################################

    pbkq = Query(exprs=dict(b=Bulks.id(), s=Bulks['stordir']()))

    bcols = ['n_atoms', 'n_elems', 'composition', 'elems',
             'lat', 'vol', 'volrat']

    abpb = PyBlock(analyze_bulks, args=[pbkq['s']], outnames=bcols)

    pop_bulks =                                                            \
        Gen(name='pop_bulks',
            desc='populates Bulks by aggregating over material+calc info',
            query=pbkq, funcs=[abpb],
            actions=[Bulks(bulks=pbkq['b'], **{x: abpb[x] for x in bcols})])

    ########################################################################
    pairs = [(Bulks[x](), Bulks['expt_'+x]()) for x in ['ce', 'bm', 'lat']]
    ss = [OR([NOT(NULL(y)), NULL(x)]) for y, x in pairs]
    ndq = Query(dict(b=Bulks.id(), s=AND(ss)),
                opt_attr=[item for sl in pairs for item in sl])
    suc = Gen(name='suc', query=ndq,
              actions=[Bulks(bulks=ndq['b'], success=ndq['s'])])
    ########################################################################
    obcols = ['bm', 'bulkmod_lstsq', 'lstsq_abc', 'eosbm']
    obq = Query(exprs=dict(
        b=Bulks.id(), s=Bulks['stordir'](), v=Bulks['volumes'](),
        e=Bulks['energies']()))
    obpb = PyBlock(analyze_opt_bulks, args=[obq[x] for x in 'sve'],
                   outnames=obcols)
    opt_bulks = Gen(
        name='opt_bulks',
        desc='populates Bulk properties from optimum bulk jobs',
        funcs=[obpb], query=obq,
        actions=[Bulks(bulks=obq['b'], **{x: obpb[x] for x in obcols})])

    ########################################################################
    eobm = Bulks['eosbm']()
    eodq = Query(dict(b=Bulks.id(),
                      d=(Bulks['bm']() - eobm)/eobm))

    eos_diff = Gen(
        name='eos_diff',
        desc='Compute Bulks.eos_diff by taking ratio of finite bulkmod '
        'to the eos bulkmod',
        query=eodq,
        actions=[Bulks(bulks=eodq['b'], eos_diff=eodq['d'])])

    ########################################################################
    cd = dict(Hydride=['H'], III=['Al'], IV=['C', 'Se', 'As'],
              V=['N', 'V', 'P'], VI=['S', 'O'], VII=['I', 'F', 'Br'])

    case = CASE([(OR([LIKE(Bulks['name'](),
                           Lit('_%{}%'.format(v))) for v in vs]
                     ), Lit(k))
                 for k, vs in cd.items()],
                Lit('Metal'))

    aq = Query(dict(b=Bulks.id(), k=case))

    alloy = Gen(
        name='alloy',
        desc='Populates bulks.alloy',
        query=aq,
        actions=[Bulks(bulks=aq['b'], alloy=aq['k'])]
    )

    ########################################################################

    iq = Query(dict(b=Bulks.id(), i=GT(ABS(Bulks['eos_diff']()), Lit(0.15))),
               opt_attr=[Bulks['expt_bm']()])
    irreg = Gen(name='irreg',
                desc='Sets Bulks.irregular',
                query=iq,
                actions=[Bulks(bulks=iq['b'], irregular=iq['i'])])
    ########################################################################
    expts = [p[1] for p in pairs]
    hq = Query(dict(b=Bulks.id(), h=NOT(AND([NULL(e) for e in expts]))),
               opt_attr=expts)
    hasdata = Gen(name='hasdata',
                  desc='True if we have any data for ce/bm/lat',
                  query=hq,
                  actions=[Bulks(bulks=hq['b'], hasdata=hq['h'])])
    ########################################################################
    ########################################################################
    ########################################################################
    gens = [pop_bulks, opt_bulks, eos_diff, alloy, irreg, struct, suc, hasdata]
    mod.add(gens)
