# External
from typing import Tuple as T
from json import loads
# Internal
from functionals.scripts.fit.get_minjob import get_minjob
from dbgen import (Model, Gen, Query, PyBlock,
                   CONCAT, MAX, Literal as Lit, GROUP_CONCAT, CASE, OR, LIKE,
                   GT, ABS)

from functionals.scripts.load.analyze_bulks        import analyze_bulks

"""
Generators purely related to analysis of Bulks and Bulkjobs
"""

def bulk_analysis(mod : Model) -> None:
    tabs = ['Bulks', 'Bulkjob', 'Calc', 'Functional']
    Bulks, Bulkjob, Calc, Fx = map(mod.get,tabs)
    ########################################################################
    ########################################################################
    ############################################################################

    biq = Query(dict(b=Bulkjob.id(), p=Bulkjob['stordir']()))

    bicols = ['parent', 'strain', 'volume', 'lattice']

    def bifunc(pth:str) -> T[str, int, float, float]:
        from os.path import join
        from ase.io import read # type: ignore
        from   numpy.linalg import norm # type: ignore
        parent, strain = pth.split('/')[-2], pth.split('/')[-1]
        atoms = read(join(pth,'POSCAR'))
        lat = norm(atoms.get_cell()[0])
        return parent, int(strain.split('_')[-1]), atoms.get_volume(), lat

    bipb = PyBlock(bifunc, args=[biq['p']], outnames=bicols)

    bulkinfo = \
        Gen(name='bulkinfo',
            desc='adds details to bulkjobs',
            funcs=[bipb],
            query=biq,
            actions=[Bulkjob(bulkjob=biq['b'], **{x:bipb[x] for x in bicols})])

    ########################################################################

    pbkqcols = ['volume','pw','fx','parent','energy','contribs','lattice', 'mag']

    vol, pw, fx, par, eng, contribs, lat, mag = [Bulkjob[x]() for x in pbkqcols]

    cats = [CONCAT([Lit('['),GROUP_CONCAT(x, order=vol),Lit(']')])
            for x in [vol,eng,contribs,lat,Bulkjob['strain'](), mag]]

    pbkq = Query(exprs = dict(n=Bulkjob['parent'](),par=par,
                              p=pw, f = fx, **dict(zip('vecltm',cats)),
                              s=MAX(Bulkjob['stordir']())),
                 opt_attr = [contribs],
                 aggcols = [par, pw, fx])

    icalc_ = Calc(insert=True, pw=pbkq['p'],
                 functional=Fx(insert=True, data=pbkq['f']))

    bcols  = ['stordir','n_atoms','n_elems','composition','elems', 'morejobs',
              'success', 'eng', 'lattice', 'volume', 'volrat', 'bulkmod',
              'bulkmod_lstsq', 'lstsq_abc', 'eosbm', 'mag',  'energies',
              'volumes',  'contribs', ]

    abpb = PyBlock(analyze_bulks,
                   args=[pbkq[x] for x in 'svelcmt'],
                   outnames=bcols)

    pop_bulks =                                                                 \
        Gen(name = 'pop_bulks',
            desc = 'populates Bulks by aggregating over material+calc info',
            query = pbkq, funcs = [abpb],
            actions = [Bulks(insert=True, name=pbkq['par'], calc=icalc_,
                             allvol=pbkq['v'], alleng=pbkq['e'],
                             strains=pbkq['t'], **{x:abpb[x] for x in bcols})])

    ########################################################################
    eobm = Bulks['eosbm']()

    eodq = Query(dict(b=Bulks.id(),
                      d=(Bulks['bulkmod']() - eobm)/eobm))

    eos_diff =                                                                  \
        Gen(name ='eos_diff',
            desc = 'Compute Bulks.eos_diff by taking ratio of finite bulkmod '
                   'to the eos bulkmod',
            query = eodq,
            actions=[Bulks(bulks=eodq['b'],eos_diff=eodq['d'])])

    ########################################################################
    cd = dict(Hydride=['H'],III=['Al'],IV=['C','Se','As'],V=['N','V','P'],
                VI=['S','O'],VII=['I','F','Br'])

    case = CASE([(OR([LIKE(Bulks['name'](),
                           Lit('_%{}%'.format(v))) for v in vs]
                   ), Lit(k))
                    for k,vs in cd.items()],
                Lit('Metal'))

    aq = Query(dict(b=Bulks.id(),k=case))

    alloy = Gen(
        name    = 'alloy',
        desc    = 'Populates bulks.alloy',
        query   = aq,
        actions = [Bulks(bulks=aq['b'],alloy=aq['k'])]
    )

    ########################################################################

    iq = Query(dict(b=Bulks.id(), i = GT(ABS(Bulks['eos_diff']()), Lit(0.15))),
            opt_attr = [Bulks['expt_bm']()])
    irreg =                                                                     \
        Gen(name = 'irreg',
            desc = 'Sets Bulks.irregular to TRUE if ASE eos bulkmod differs from discrete by too much',
            query= iq,
            actions=[Bulks(bulks=iq['b'],irregular=iq['i'])])
    ########################################################################


    ########################################################################
    ########################################################################
    ########################################################################
    gens = [bulkinfo, pop_bulks, eos_diff, alloy, irreg]
    mod.add(gens)
