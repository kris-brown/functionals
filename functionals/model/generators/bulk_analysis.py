# Internal
from dbgen import (Model, Gen, Query, PyBlock, NOT, AND,
                   Literal as Lit, CASE, OR, LIKE, NULL,
                   EQ)

from functionals.scripts.load.analyze_bulks import analyze_bulks
from functionals.scripts.load.analyze_opt_bulks import analyze_opt_bulks
from functionals.scripts.load.primcell import primcell


"""
Generators purely related to analysis of Bulks and Bulkjobs.
"""


def bulk_analysis(mod: Model) -> None:
    tabs = ['Bulks', 'Calc']
    Bulks, Calc = map(mod.get, tabs)
    ########################################################################
    ########################################################################
    ########################################################################

    pbkq = Query(exprs=dict(b=Bulks.id(), s=Bulks['stordir'](),
                            t=Bulks['name']()))

    bcols = ['n_atoms', 'n_elems', 'composition', 'elems', 'struct']

    abpb = PyBlock(analyze_bulks, args=[pbkq['s'], pbkq['t']], outnames=bcols)

    pop_bulks =                                                            \
        Gen(name='pop_bulks',
            desc='populates Bulks by aggregating over material+calc info',
            query=pbkq, transforms=[abpb],
            loads=[Bulks(bulks=pbkq['b'], **{x: abpb[x] for x in bcols})])

    ########################################################################
    pairs = [(Bulks[x](), Bulks['expt_' + x]()) for x in ['ce', 'bm', 'lat']]
    ss = [OR([NOT(NULL(y)), NULL(x)]) for y, x in pairs]
    ndq = Query(dict(b=Bulks.id(), s=AND(ss)),
                opt_attr=[item for sl in pairs for item in sl])
    suc = Gen(name='suc', query=ndq,
              loads=[Bulks(bulks=ndq['b'], success=ndq['s'])])
    ########################################################################
    obcols = ['bm', 'cellvol', 'bulkmod_lstsq', 'lstsq_abc']
    obq = Query(exprs=dict(
        b=Bulks.id(), s=Bulks['stordir'](), v=Bulks['volumes'](),
        e=Bulks['energies'](), t=Bulks['struct']()),
        constr=EQ(Bulks['err'](), Lit('')))
    obpb = PyBlock(analyze_opt_bulks, args=[obq[x] for x in 'svet'],
                   outnames=obcols)
    opt_bulks = Gen(
        name='opt_bulks',
        desc='populates Bulk properties from optimum bulk jobs',
        transforms=[obpb], query=obq,
        loads=[Bulks(bulks=obq['b'], **{x: obpb[x] for x in obcols})])

    ########################################################################
    vvq = Query(dict(b=Bulks.id(), p=Bulks['stordir'](),
                     v=Bulks['vol'](), r=Bulks['volrat'](),
                     s=Bulks['struct']()))

    def get_lat(st: str, v: float, vr: float, pth: str) -> float:
        a = (float(v) / float(vr))**(1 / 3)
        if 'prim' not in st:
            return a
        elif True:  # 'zincblende' in st or 'rocksalt' in st:
            return a * 2 ** 0.5
        else:
            breakpoint()
            raise ValueError()

    glpb = PyBlock(get_lat, args=[vvq[x] for x in 'svrp'])
    compute_lat = Gen(
        name='compute_lat',
        desc='Compute DFT lattice constant from optimum volume and vol ratio',
        query=vvq, transforms=[glpb],
        loads=[Bulks(bulks=vvq['b'], lat=glpb['out'])])

    ########################################################################
    pcq = Query(dict(b=Bulks.id(), p=Bulks['stordir'](), s=Bulks['struct'](),
                     v=Bulks['cellvol'](), r=Bulks['volprimrat']()))
    pccols = ['volrat', 'vol', 'prim']
    pcpb = PyBlock(primcell, args=[pcq[x] for x in 'psvr'], outnames=pccols)
    primcell_ = Gen(
        name='primcell',
        desc='Compute DFT lattice constant from optimum volume and vol ratio',
        query=pcq, transforms=[pcpb],
        loads=[Bulks(bulks=pcq['b'], **{x: pcpb[x] for x in pccols})])
    ########################################################################
    vpq = Query(dict(b=Bulks.id(), s=Bulks['struct']()))
    vppb = PyBlock(
        (lambda z: {**{"": 0}, **dict(
            fcc=4, bcc=2, rocksalt=4, diamond=4, zincblende=4, hcp=1,
            cesiumchloride=1, wurtzite=0)}[z]),
        args=[vpq['s']])
    vprat = Gen(
        name='volprimrat',
        desc='Look up ratio of primative volume to conventional volume',
        query=vpq, transforms=[vppb],
        loads=[Bulks(bulks=vpq['b'], volprimrat=vppb['out'])])
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
        loads=[Bulks(bulks=aq['b'], alloy=aq['k'])]
    )

    ########################################################################
    expts = [p[1] for p in pairs]
    hq = Query(dict(b=Bulks.id(), h=NOT(AND([NULL(e) for e in expts]))),
               opt_attr=expts)
    hasdata = Gen(name='hasdata',
                  desc='True if we have any data for ce/bm/lat',
                  query=hq,
                  loads=[Bulks(bulks=hq['b'], hasdata=hq['h'])])
    ########################################################################
    # gens = []  # type: List[Gen]
    # for k in ['ce', 'bm', 'lat']:
    #     eq = Query(dict(b=Bulks.id(), e=Bulks[k]()-Bulks['expt_'+k]()))
    #     gens.append(Gen(name='err '+k, query=eq, loads=[
    #         Bulks(bulks=eq['b'], **{'err_'+k: eq['e']})]))
    # ec, eb, el = gens
    ########################################################################
    ########################################################################
    ########################################################################
    gens = [pop_bulks, opt_bulks, alloy, suc, hasdata, primcell_, vprat,
            compute_lat]

    mod.add(gens)
