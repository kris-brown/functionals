# External
from typing import Tuple as T
from json import loads, dumps
from ast import literal_eval
import numpy as np

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, Import,
                   Literal as Lit, EQ, LEFT, One, SUM,
                   CONCAT, AVG, LIKE, GROUP_CONCAT, ABS)

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
               query=tmagq, funcs=[tmagpb],
               actions=[Atoms(atoms=tmagq['i'], true_mag=tmagpb['out'])])

    ########################################################################
    xvq = Query(dict(b=Bulks.id(),
                     v=(Bulks['expt_lat']()**Lit(3) * Bulks['volrat']())))
    xvol = Gen('xvol', query=xvq, actions=[Bulks(bulks=xvq['b'],
                                                 expt_vol=xvq['v'])])
    ########################################################################
    ibq = Query(exprs=dict(c=Calc.id(),
                           b=EQ(LEFT(Calc['data'](), One), Lit('['))))
    isbeef =                                                                  \
        Gen(name='isbeef',
            desc='Check whether fx was BEEF - if fx is a word or array',
            query=ibq,
            actions=[Calc(calc=ibq['c'], beef=ibq['b'])])

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
            query=erq, funcs=[erpb],
            actions=[Refs(insert=True,
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
            query=eq, actions=[Bulks(bulks=eq['b'], eform=eq['e'])])

    ########################################################################

    ceq = Query(exprs=dict(b=Bulks.id(),
                           ce=Lit(-1) * Bulks['eform']() / Bulks['n_atoms']()))
    ce =                                                                     \
        Gen(name='ce',
            desc='Difference between a relaxed energy (per atom) and '
                 'the sum of reference energies (divided by # of atoms)',
            query=ceq,
            actions=[Bulks(bulks=ceq['b'], ce=ceq['ce'])])

    ########################################################################

    bpth = mod.make_path('bulks', [bulks__calc])
    com = Lit(',')
    bname, bce, bbm = [Bulks[x](bpth) for x in ['name', 'ce', 'bm']]
    catargs = [bname, com, bce, com, bbm]
    dq = Query(exprs=dict(c=Calc.id(), w=Calc['allmat'](),
                          h=GROUP_CONCAT(CONCAT(catargs), delim=' ')),
               basis=[Calc], aggcols=[Calc.id()], opt_attr=[bce, bbm])

    def dfun(have_: str, want_: str) -> T[str, int, bool]:
        wce, wbml = [set(x.split()) for x in want_.split('|')]

        for name, ce, bm in [x.split(',') for x in have_.split()]:
            if ce and name in wce:
                wce.remove(name)
            if bm and name in wbml:
                wbml.remove(name)
        missing = [x + '_ce' for x in wce] + [x + '_bml' for x in wbml]
        n = len(missing)
        return ' '.join(sorted(missing)), n, n == 0

    dcols = ['missing', 'n_missing', 'done']
    dqpb = PyBlock(dfun, args=[dq['h'], dq['w']], outnames=dcols)

    done =                                                                  \
        Gen(name='done',
            desc='Determines if calculator is ready for fitting (has data)',
            query=dq, funcs=[dqpb],
            actions=[Calc(calc=dq['c'], **{x: dqpb[x] for x in dcols})])

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

    refcontribs =                                                   \
        Gen(name='refcontribs',
            desc='Takes atomic contribs and multiplies by stoichiometry',
            query=rcq,
            funcs=[rcpb],
            actions=[Refs(refs=rcq['r'],
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
               actions=[Bulks(**{x: errq[x] for x in errq.exprs.keys()})])
    ########################################################################
    # cols = dict(ce=('ce', 'expt_ce'), bm=('bm', 'expt_bm'),
    #             lat=('lat', 'expt_lat'),  mag=('mag', 'expt_mag'))
    cbp = mod.make_path('bulks', [bulks__calc])
    msegens = []
    for m in mets:
        mae, relmae = [Bulks.get(x + 'err_' + m)(cbp) for x in ['', 'rel']]

        mseq = Query(exprs=dict(c=Calc.id(), m=AVG(ABS(mae)),
                                r=AVG(ABS(relmae))),
                     basis=[Calc],
                     aggcols=[Calc.id()])
        msegens.append(Gen(
            name='mae_' + m, query=mseq, tags=['mae'],
            actions=[Calc(calc=mseq['c'], **{
                'mae_' + m: mseq['m'], 'relmae_' + m:mseq['r']})]))
    msec, mseb, msel, msev, msem = msegens
    #####################################################################

    annq = Query(dict(a=Atoms.id(), s=Atoms['stordir']()))

    def annfunc(pth: str) -> T[str, int]:
        from ase.data import chemical_symbols
        symb = pth.split('/')[-1]
        return symb, chemical_symbols.index(symb)

    annpb = PyBlock(annfunc, args=[annq['s']], outnames=['n', 'm'])
    atom_namenum =                                                       \
        Gen(name='atom_namenum',
            desc='populate atoms.name and atoms.num',
            query=annq,
            funcs=[annpb],
            actions=[Atoms(atoms=annq['a'], name=annpb['n'], num=annpb['m'])])

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
    #         import pdb
    #         pdb.set_trace()
    #         assert False
    #         return ''
    #     else:
    #         return data

    # cnpb = PyBlock(cnfunc, args=[cnq['d']])
    # calcname2 =                                                           \
    #     Gen(name='calcname2',
    #         desc='Figures out the name of a BEEF calculator from its coefs',
    #         query=cnq, funcs=[cnpb],
    #         actions=[Calc(calc=cnq['c'], name=cnpb['out'])])
    ########################################################################
    seq = Query(dict(s=Surf.id(), m=Surf['mat'](
    ), o=Surf['ontop'](), h=Surf['hollow']()))

    def sef(m: str, o: float, h: float) -> float:
        if m == 'Pd':
            return h - o
        else:
            return o - h
    sepb = PyBlock(sef, args=[seq[x] for x in 'moh'])
    surfan = \
        Gen(name='surferr',
            query=seq, funcs=[sepb], tags=['surf'],
            actions=[Surf(surf=seq['s'], err=sepb['out'])])
    ########################################################################
    ########################################################################
    ########################################################################

    gens = [isbeef, exptref, eform, ce, refcontribs, xvol, msem,
            msec, mseb, msel, msev, atom_namenum,
            tmag, done, surfan, errs]

    mod.add(gens)
