# External
from typing import Tuple as T
# Internal
from dbgen import (Model, Gen, Query, PyBlock, JPath,
                   CONCAT, MAX, Literal as Lit, GROUP_CONCAT, CASE, OR, LIKE,
                   GT, ABS, NOT, AND, NULL, LEFT, POSITION, EQ)

from functionals.scripts.load.analyze_bulks import analyze_bulks
from functionals.scripts.load.analyze_opt_bulks import analyze_opt_bulks


"""
Generators purely related to analysis of Bulks and Bulkjobs.
"""


def bulk_analysis(mod: Model) -> None:
    tabs = ['Bulks', 'Bulkjob', 'Calc', 'Functional']
    Bulks, Bulkjob, Calc, Fx = map(mod.get, tabs)
    bulks__calc = mod.get_rel(Bulks.r('calc'))
    calc__functional = mod.get_rel(Calc.r('functional'))
    ########################################################################
    ########################################################################
    ########################################################################

    def get_struct(mat: str) -> str:
        from functionals.CLI.submit import matdata
        try:
            return str(matdata[mat].struct)
        except KeyError:
            print('\n\n\n\n', mat)
            import pdb
            pdb.set_trace()
            return ''

    sq = Query(dict(b=Bulks.id(), n=Bulks['name']()))
    sf = PyBlock(get_struct, args=[sq['n']])
    struct = \
        Gen(name='struct', query=sq, funcs=[sf],
            actions=[Bulks(bulks=sq['b'], struct=sf['out'])])
    ########################################################################

    biq = Query(dict(b=Bulkjob.id(), p=Bulkjob['stordir']()))

    bcols = ['parent', 'strain', 'volume', 'lattice', 'opt']

    def bifunc(pth: str) -> T[str, int, float, float, bool]:
        from os.path import join
        from ase.io import read
        from numpy.linalg import norm
        parent, strain = pth.split('/')[-2], pth.split('/')[-1]
        atoms = read(join(pth, 'POSCAR'))
        lat = norm(atoms.get_cell()[0])
        s, o = int(strain.split('_')[-1]), 'optstrain' in pth
        return parent, s, atoms.get_volume(), lat, o

    bipb = PyBlock(bifunc, args=[biq['p']], outnames=bcols)

    bulkinfo = \
        Gen(name='bulkinfo',
            desc='adds details to bulkjobs',
            funcs=[bipb],
            query=biq,
            actions=[Bulkjob(bulkjob=biq['b'], **{x: bipb[x] for x in bcols})])

    ########################################################################

    pbkqcol = ['volume', 'pw', 'fx', 'parent', 'energy', 'contribs',
               'lattice', 'mag']

    vol, pw, fx, par, eng, contribs, lat, mag = [Bulkjob[x]() for x in pbkqcol]

    cats = [CONCAT([Lit('['), GROUP_CONCAT(x, order=vol), Lit(']')])
            for x in [vol, eng, Bulkjob['strain'](), mag, contribs, lat]]
    mxs = MAX(Bulkjob['stordir']())
    pbkq = Query(exprs=dict(n=Bulkjob['parent'](), par=par,
                            p=pw, f=fx, **dict(zip('vetmcl', cats)),
                            s=mxs),
                 constr=NOT(Bulkjob['opt']()),
                 opt_attr=[contribs],
                 aggcols=[par, pw, fx])

    icalc_ = Calc(insert=True, pw=pbkq['p'],
                  functional=Fx(insert=True, data=pbkq['f']))

    bcols = ['stordir', 'n_atoms', 'n_elems', 'composition', 'elems',
             'morejobs', 'success', 'eosbm', 'lat', 'vol', 'volrat']

    abpb = PyBlock(analyze_bulks,
                   args=[pbkq[x] for x in 'svelcmt'],
                   outnames=bcols)

    pop_bulks =                                                            \
        Gen(name='pop_bulks',
            desc='populates Bulks by aggregating over material+calc info',
            query=pbkq, funcs=[abpb],
            actions=[Bulks(insert=True, name=pbkq['par'], calc=icalc_,
                           allvol=pbkq['v'], alleng=pbkq['e'],
                           strains=pbkq['t'], **{x: abpb[x] for x in bcols})])

    ########################################################################

    mets = [Bulks['expt_'+x]() for x in ['ce', 'bm', 'lat']]
    ndq = Query(dict(b=Bulks.id(), r=AND([NULL(x) for x in mets])),
                opt_attr=mets)
    nd = Gen(name='nodata', query=ndq,
             actions=[Bulks(bulks=ndq['b'], nodata=ndq['r'])])
    ########################################################################
    obcols = ['bm', 'bulkmod_lstsq', 'lstsq_abc', 'mag',
              'energy', 'optsuccess']
    # cats = vol, eng, strain, mag, contribs, lat
    extra = dict(zip('vetmclxyz', cats+[NOT(NULL(MAX(x))) for x in mets]))
    Len = POSITION(Lit('optstrain'), mxs) - Lit(2)
    trunc = LEFT(mxs, Len)
    jp = JPath(Fx, [calc__functional, bulks__calc])
    obq = Query(exprs=dict(n=Bulkjob['parent'](), par=par,
                           p=pw, f=fx, **extra, s=trunc),
                basis=[Bulkjob],
                constr=AND([Bulkjob['opt'](),
                            EQ(Fx['data'](jp), Bulkjob['fx'](),),
                            EQ(Bulks['name'](), Bulkjob['parent']())]),
                opt_attr=[contribs]+mets,
                aggcols=[par, pw, fx])

    icalc_ = Calc(pw=obq['p'],
                  functional=Fx(data=obq['f']))

    obpb = PyBlock(analyze_opt_bulks, args=[obq[x] for x in 'vemtxyzn'],
                   outnames=obcols)
    opt_bulks =                                                            \
        Gen(name='opt_bulks',
            desc='populates Bulk properties from optimum bulk jobs',
            funcs=[obpb], query=obq,
            actions=[Bulks(stordir=obq['s'], calc=icalc_,
                           volumes=obq['v'], energies=obq['e'],
                           contribs=obq['c'], optstrains=obq['t'],
                           **{x: obpb[x] for x in obcols})])

    ########################################################################
    eobm = Bulks['eosbm']()
    eodq = Query(dict(b=Bulks.id(),
                      d=(Bulks['bm']() - eobm)/eobm))

    eos_diff =                                                             \
        Gen(name='eos_diff',
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
    irreg =                                                                \
        Gen(name='irreg',
            desc='Sets Bulks.irregular',
            query=iq,
            actions=[Bulks(bulks=iq['b'], irregular=iq['i'])])
    ########################################################################
    ########################################################################
    ########################################################################
    gens = [bulkinfo, pop_bulks, opt_bulks, eos_diff, alloy, irreg, nd, struct]
    mod.add(gens)
