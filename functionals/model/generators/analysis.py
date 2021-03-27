# External
from typing import Tuple as T
from json import loads, dumps
from ast import literal_eval
import numpy as np

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, Import,
                   Literal as Lit, EQ, LEFT, One, SUM,
                   CONCAT, AVG, LIKE, ABS, COUNT)

from functionals.scripts.load.true_mag import true_mag

#########################################################################
#########################################################################
#########################################################################


def analysis(mod: Model) -> None:
    # Extract tables
    tabs = ['Bulks', 'Atoms', 'Refs', 'Calc', 'Surf']
    Bulks, Atoms, Refs, Calc, Surf = map(mod.get, tabs)

    refs__atoms, refs__bulks, bulks__calc, atoms__calc = \
        map(mod.get_rel, [Refs.r('atoms'), Refs.r('bulks'), Bulks.r('calc'),
                          Atoms.r('calc')])
    ########################################################################
    ########################################################################

    ########################################################################
    tmagq = Query(dict(i=Atoms.id(), n=Atoms['name']()))
    tmagpb = PyBlock(true_mag, args=[tmagq['n']])
    tmag = Gen(name='truemag',
               desc='populates atoms.true_mag',
               query=tmagq, transforms=[tmagpb],
               loads=[Atoms(atoms=tmagq['i'], true_mag=tmagpb['out'])])

    ########################################################################
    xvq = Query(dict(b=Bulks.id(),
                     v=(Bulks['expt_lat']()**Lit(3) * Bulks['volrat']())))
    xvol = Gen('xvol', query=xvq, loads=[Bulks(bulks=xvq['b'],
                                               expt_vol=xvq['v'])])
    ########################################################################
    ibq = Query(exprs=dict(c=Calc.id(),
                           b=EQ(LEFT(Calc['data'](), One), Lit('['))))
    isbeef =                                                                  \
        Gen(name='isbeef',
            desc='Check whether fx was BEEF - if fx is a word or array',
            query=ibq,
            loads=[Calc(calc=ibq['c'], beef=ibq['b'])])

    #####################################################################
    name = CONCAT([Lit('%,'), Atoms['num'](), Lit(',%')])

    bcalc, acalc = [mod.make_path('calc', [x]) for x in
                    [bulks__calc, atoms__calc]]

    erq = Query(exprs=dict(b=Bulks.id(), a=Atoms.id(),
                           c=Bulks['composition'](),
                           n=Atoms['num']()),
                basis=[Bulks, Atoms],
                constr=AND([
                    EQ(Calc.id(bcalc), Calc.id(acalc)),
                    LIKE(Bulks['elems'](), name),
                    EQ(Atoms['err'](), Lit(''))]))

    erpb = PyBlock(lambda c, n: literal_eval(c).get(n, []),
                   env=Env([Import('ast', ['literal_eval'])]),
                   args=[erq['c'], erq['n']])

    exptref =                                                                \
        Gen(name='refs',
            desc='Populates expt_refs mapping table',
            query=erq, transforms=[erpb],
            loads=[Refs(insert=True,
                        bulks=erq['b'],
                        atoms=erq['a'],
                        num=erpb['out'])])

    ########################################################################
    # Paths
    refpth = mod.make_path("refs", [refs__bulks])

    eq = Query(exprs=dict(b=Bulks.id(),
                          e=Bulks['eng']() - SUM(Refs['energy'](refpth))),
               aggcols=[Bulks.id()],
               basis=[Bulks],
               aconstr=EQ(SUM(Refs['num'](refpth)), Bulks['n_atoms']()))

    eform =                                                                  \
        Gen(name='eform',
            desc='Diff between a relaxed energy the sum of reference engs ',
            query=eq, loads=[Bulks(bulks=eq['b'], eform=eq['e'])])

    ########################################################################

    ceq = Query(exprs=dict(b=Bulks.id(),
                           ce=Lit(-1) * Bulks['eform']() / Bulks['n_atoms']()))
    ce =                                                                     \
        Gen(name='ce',
            desc='Difference between a relaxed energy (per atom) and '
                 'the sum of reference energies (divided by # of atoms)',
            query=ceq,
            loads=[Bulks(bulks=ceq['b'], ce=ceq['ce'])])

    ########################################################################

    bpth = mod.make_path('bulks', [bulks__calc])
    bname, bce, bbm, bsucc = [Bulks[x](bpth)
                              for x in ['name', 'ce', 'bm', 'success']]
    n_present = COUNT(Lit(1))
    l116 = Lit(116)
    dq = Query(exprs=dict(c=Calc.id(),
                          n=l116 - n_present,
                          d=EQ(n_present, l116)),
               basis=[Calc], aggcols=[Calc.id()], constr=bsucc)

    done =                                                                  \
        Gen(name='done',
            desc='Determines if calculator is ready for fitting (has data)',
            query=dq,
            loads=[Calc(calc=dq['c'], done=dq['d'], n_missing=dq['n'])])

    ########################################################################
    refatom = mod.make_path('atoms', [refs__atoms])

    rcq = Query(exprs=dict(r=Refs.id(),
                           c=Atoms['contribs'](refatom),
                           n=Refs['num'](),
                           e=Refs['num']() * Atoms['energy'](refatom)),
                basis=[Refs],
                opt_attr=[Atoms['contribs'](refatom)])

    rcpb = PyBlock(
        lambda c, n: dumps((n * np.array(loads(c))).tolist()) if c else None,
        args=[rcq[x] for x in 'cn'])

    refcontribs = Gen(name='refcontribs',
                      desc='Takes atomic contribs and multiplies by stoich',
                      query=rcq,
                      transforms=[rcpb],
                      loads=[Refs(refs=rcq['r'],
                                  energy=rcq['e'],
                                  contribs=rcpb['out'])])

    ########################################################################
    mets = ['ce', 'bm', 'lat', 'vol', 'mag']
    errcols = {x + 'err_' + k: (Bulks.get(k)() - Bulks.get('expt_' + k)()) / (
        Bulks.get('expt_' + k)() if x else Lit(1)) for x in ['rel', '']
        for k in mets}
    errq = Query(dict(bulks=Bulks.id(), **errcols),
                 opt_attr=[Bulks.get(x + k)() for x in ['', 'expt_']
                           for k in mets])

    errs = Gen('errs', desc='Difference btw computed and expt',
               query=errq,
               loads=[Bulks(**{x: errq[x] for x in errq.exprs.keys()})])
    ########################################################################
    # cols = dict(ce=('ce', 'expt_ce'), bm=('bm', 'expt_bm'),
    #             lat=('lat', 'expt_lat'),  mag=('mag', 'expt_mag'))
    cbp = mod.make_path('bulks', [bulks__calc])
    msegens = []
    for m in mets:
        mae, relmae = [Bulks.get(x + 'err_' + m)(cbp) for x in ['', 'rel']]

        mseq = Query(exprs=dict(c=Calc.id(), m=AVG(ABS(mae)),
                                r=AVG(ABS(relmae))),
                     constr=EQ(Bulks['err'](cbp), Lit("")),
                     basis=[Calc],
                     aggcols=[Calc.id()])
        msegens.append(Gen(
            name='mae_' + m, query=mseq, tags=['mae'],
            loads=[Calc(calc=mseq['c'], **{
                'mae_' + m: mseq['m'], 'relmae_' + m:mseq['r']})]))
    msec, mseb, msel, msev, msem = msegens
    #####################################################################

    annq = Query(dict(a=Atoms.id(), s=Atoms['stordir']()))

    def annfunc(pth: str) -> T[str, int]:
        from ase.data import chemical_symbols
        symb = pth.split('/')[-1]
        return symb, chemical_symbols.index(symb)

    annpb = PyBlock(annfunc, args=[annq['s']], outnames=['n', 'm'])
    atom_namenum = Gen(name='atom_namenum',
                       desc='populate atoms.name and atoms.num',
                       query=annq,
                       transforms=[annpb],
                       loads=[Atoms(atoms=annq['a'], name=annpb['n'],
                                    num=annpb['m'])])

    ########################################################################

    # cnq = Query(dict(c=Calc.id(), d=Fx['data'](cnp)), basis=[Calc])

    # def cnfunc(data: str) -> str:
    #     import os
    #     import json
    #     root = '/Users/ksb/functionals/data/beefs'
    #     if data[0] == '[':
    #         data_ = json.loads(data)
    #         for beef in os.listdir(root):
    #             with open(os.path.join(root, beef), 'r') as f:
    #                 if json.load(f) == data_:
    #                     return beef
    #
    #         breakpoint()
    #         assert False
    #         return ''
    #     else:
    #         return data

    # cnpb = PyBlock(cnfunc, args=[cnq['d']])
    # calcname2 =                                                           \
    #     Gen(name='calcname2',
    #         desc='Figures out the name of a BEEF calculator from its coefs',
    #         query=cnq, transforms=[cnpb],
    #         loads=[Calc(calc=cnq['c'], name=cnpb['out'])])
    ########################################################################
    seq = Query(dict(s=Surf.id(), m=Surf['mat'](
    ), o=Surf['ontop'](), h=Surf['hollow']()))

    def sef(m: str, o: float, h: float) -> float:
        if m == 'Pd':
            return h - o
        else:
            return o - h
    sepb = PyBlock(sef, args=[seq[x] for x in 'moh'])
    surfan = Gen(name='surferr',
                 query=seq, transforms=[sepb], tags=['surf'],
                 loads=[Surf(surf=seq['s'], err=sepb['out'])])
    ########################################################################
    ########################################################################
    ########################################################################

    gens = [isbeef, exptref, eform, ce, refcontribs, xvol, msem,
            msec, mseb, msel, msev, atom_namenum,
            tmag, done, surfan, errs]

    mod.add(gens)
